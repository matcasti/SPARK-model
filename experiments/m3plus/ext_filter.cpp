// EXPERIMENTAL generalized filter for M3 and extensions (not part of the pipeline; verified against engine/fast_filter.cpp in lib.R).
// State S_k = [x_p, x_s, z1, z2, I_net]: the two OU drives, an optional damped oscillator (respiratory/baroreflex) and the net phase integral
// I_net = -I_p + I_s - I_o accumulated over the interval that ends at beat k. The transition has a zero column for I_net (the reset), so the
// recursion is an ordinary Kalman filter and RTS smoothing applies directly.
// Streams: stream 0 has target 1 (phase); streams 1..m are regressors (column 0 of Xreg = dt gives nu0), profiled by GLS.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <complex>
#include <Eigen/SVD>
using namespace Rcpp;
typedef Eigen::Matrix<double, 5, 5> M5; typedef Eigen::Matrix<double, 5, 1> V5;
static M5 pinv5(const M5& X) {
  Eigen::JacobiSVD<M5> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
  V5 sv = svd.singularValues(); double tol = 1e-10 * std::max(sv(0), 1e-300); M5 Si = M5::Zero();
  for (int i = 0; i < 5; ++i) if (sv(i) > tol) Si(i, i) = 1.0 / sv(i);
  return svd.matrixV() * Si * svd.matrixU().transpose();
}

struct OU { double phi, q11, q13, q33; };
static inline OU ou(double k, double dt, double V) {            // exact OU integrals (mirrors .compute_ou_block)
  double x = k * dt, s2 = 2.0 * k * V; OU b;
  if (x < 1e-3) {
    b.phi = dt - (k * dt * dt) / 2.0 + (k * k * dt * dt * dt) / 6.0;
    b.q11 = s2 * (dt - k * dt * dt + (2.0 * k * k * std::pow(dt, 3)) / 3.0);
    b.q13 = s2 * ((dt * dt) / 2.0 - (k * std::pow(dt, 3)) / 2.0 + (7.0 * k * k * std::pow(dt, 4)) / 24.0);
    b.q33 = s2 * (std::pow(dt, 3) / 3.0 - (k * std::pow(dt, 4)) / 4.0 + (7.0 * k * k * std::pow(dt, 5)) / 60.0);
  } else {
    double e1 = std::exp(-x), e2 = std::exp(-2.0 * x);
    b.phi = (1.0 - e1) / k; b.q11 = s2 * (1.0 - e2) / (2.0 * k);
    b.q13 = s2 * (1.0 - 2.0 * e1 + e2) / (2.0 * k * k); b.q33 = s2 * (2.0 * x - 3.0 + 4.0 * e1 - e2) / (2.0 * k * k * k);
  }
  return b;
}
// Damped oscillator z = z1 + i z2, dz = (-alpha + i omega) z dt + sigma dW_c (circular complex noise), stationary Var(z1) = Var(z2) = Vo
// (sigma^2 = 4 alpha Vo). Over h: z' = e^{lambda h} z + w1, integral of z1 = Re(g z) + Re(w2), g = (e^{lambda h} - 1) / lambda.
struct Osc { double ea, c, s, gre, gim, A, Bre, Bim, C; };
static inline Osc osc_block(double alpha, double omega, double h, double Vo) {
  typedef std::complex<double> cd; cd lam(-alpha, omega), e = std::exp(lam * h), g = (e - 1.0) / lam;
  double base = (1.0 - std::exp(-2.0 * alpha * h)) / (2.0 * alpha), s2 = 4.0 * alpha * Vo;
  cd B = s2 / std::conj(lam) * (base - g);
  Osc o; o.ea = std::exp(-alpha * h); o.c = std::cos(omega * h); o.s = std::sin(omega * h); o.gre = g.real(); o.gim = g.imag();
  o.A = s2 * base; o.Bre = B.real(); o.Bim = B.imag(); o.C = s2 / std::norm(lam) * (base - 2.0 * g.real() + h);
  return o;
}
//' Transition and process covariance of the oscillator block over h (for tests): rows/cols z1, z2, I_o.
// [[Rcpp::export]]
List osc_moments(double alpha, double omega, double h, double Vo) {
  Osc o = osc_block(alpha, omega, h, Vo);
  NumericMatrix F(3, 3), Q(3, 3);
  F(0, 0) = o.ea * o.c; F(0, 1) = -o.ea * o.s; F(1, 0) = o.ea * o.s; F(1, 1) = o.ea * o.c; F(2, 0) = o.gre; F(2, 1) = -o.gim;
  Q(0, 0) = Q(1, 1) = o.A / 2; Q(0, 2) = Q(2, 0) = o.Bre / 2; Q(1, 2) = Q(2, 1) = o.Bim / 2; Q(2, 2) = o.C / 2;
  return List::create(_["F"] = F, _["Q"] = Q);
}

// [[Rcpp::export]]
List ext_filter(NumericVector dy_, NumericVector gm_, NumericVector qs_, NumericMatrix Xreg, double kp, double ks, double lp, double lR,
                double osc_alpha, double osc_omega, double osc_lo, double lik_gate_max, bool full, bool smooth,
                NumericVector beta_fixed, double sig2_fixed) {
  const int N = dy_.size(), m = Xreg.ncol(), Mst = m + 1;
  const bool use_osc = osc_lo > 0 && osc_alpha > 0;
  std::vector<V5> Xf(Mst, V5::Zero());
  M5 P = M5::Zero(); P(0, 0) = lp; P(1, 1) = 1.0; if (use_osc) { P(2, 2) = osc_lo; P(3, 3) = osc_lo; }
  std::vector<double> S(N), Sung(N); std::vector<std::vector<double> > v(Mst, std::vector<double>(N));
  bool keep = full || smooth;
  std::vector<std::vector<V5> > xfilt, xpred; std::vector<M5> Pfilt, Ppred, Fk;
  if (keep) { xfilt.assign(Mst, std::vector<V5>(N)); xpred.assign(Mst, std::vector<V5>(N)); Pfilt.resize(N); Ppred.resize(N); Fk.resize(N); }
  M5 I5 = M5::Identity();
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], qs = qs_[k], R = std::min(lR * gm_[k], 1e12) * dt;
    OU bp = ou(kp, dt, lp), bs = ou(ks, dt, 1.0);
    M5 F = M5::Zero(); M5 Q = M5::Zero();
    F(0, 0) = std::exp(-kp * dt); F(1, 1) = std::exp(-ks * dt); F(4, 0) = -bp.phi; F(4, 1) = bs.phi;
    Q(0, 0) = bp.q11 * qs; Q(0, 4) = Q(4, 0) = -bp.q13 * qs; Q(1, 1) = bs.q11 * qs; Q(1, 4) = Q(4, 1) = bs.q13 * qs; Q(4, 4) = (bp.q33 + bs.q33) * qs;
    if (use_osc) {
      Osc o = osc_block(osc_alpha, osc_omega, dt, osc_lo);
      F(2, 2) = o.ea * o.c; F(2, 3) = -o.ea * o.s; F(3, 2) = o.ea * o.s; F(3, 3) = o.ea * o.c; F(4, 2) = -o.gre; F(4, 3) = o.gim;
      Q(2, 2) = o.A / 2 * qs; Q(3, 3) = o.A / 2 * qs; Q(2, 4) = Q(4, 2) = -o.Bre / 2 * qs; Q(3, 4) = Q(4, 3) = -o.Bim / 2 * qs; Q(4, 4) += o.C / 2 * qs;
    }
    M5 Pp = F * P * F.transpose() + Q;
    double Sk = std::max(Pp(4, 4) + R, 1e-12); S[k] = Sk; Sung[k] = std::max(Pp(4, 4) + lR * dt, 1e-12);
    V5 K = Pp.col(4) / Sk;
    for (int j = 0; j < Mst; ++j) {
      V5 xp = F * Xf[j]; double target = (j == 0) ? 1.0 : Xreg(k, j - 1);
      double vj = target - xp(4); v[j][k] = vj; if (keep) xpred[j][k] = xp;
      Xf[j] = xp + K * vj; if (keep) xfilt[j][k] = Xf[j];
    }
    M5 IKH = I5; IKH.col(4) -= K;                                           // I - K H, H = e5'
    P = IKH * Pp * IKH.transpose() + K * K.transpose() * R; P = 0.5 * (P + P.transpose());
    if (keep) { Pfilt[k] = P; Ppred[k] = Pp; Fk[k] = F; }
  }
  // GLS profile of the regression coefficients (stream 0 is the response, streams 1..m the regressors)
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); Eigen::VectorXd b = Eigen::VectorXd::Zero(m);
  for (int k = 0; k < N; ++k) for (int i = 0; i < m; ++i) { b(i) += v[i + 1][k] * v[0][k] / S[k]; for (int j = 0; j < m; ++j) A(i, j) += v[i + 1][k] * v[j + 1][k] / S[k]; }
  Eigen::VectorXd beta(m);
  if (beta_fixed.size() == m) { for (int i = 0; i < m; ++i) beta(i) = beta_fixed[i]; }
  else {
    bool fin = A.allFinite() && b.allFinite();
    if (fin) beta = A.colPivHouseholderQr().solve(b);
    if (!fin || !beta.allFinite() || !(beta(0) > 0)) return List::create(_["LL_star"] = R_NegInf, _["nu0"] = NA_REAL, _["sig2"] = NA_REAL);
  }
  std::vector<double> vk(N); int n_eff = 0; for (int k = 0; k < N; ++k) if (gm_[k] < lik_gate_max) ++n_eff;
  bool all = n_eff < 2; if (all) n_eff = N;
  double s2 = 0, sumlogS = 0;
  for (int k = 0; k < N; ++k) { double r = v[0][k]; for (int i = 0; i < m; ++i) r -= beta(i) * v[i + 1][k]; vk[k] = r;
    if (all || gm_[k] < lik_gate_max) { s2 += r * r / S[k]; sumlogS += std::log(S[k]); } }
  double sig2 = (sig2_fixed > 0) ? sig2_fixed : std::max(s2 / n_eff, 1e-12);
  double LL = -0.5 * sumlogS - (n_eff / 2.0) * std::log(sig2);
  NumericVector betaR(m); for (int i = 0; i < m; ++i) betaR[i] = beta(i);
  List out = List::create(_["LL_star"] = LL, _["nu0"] = beta(0), _["sig2"] = sig2, _["beta"] = betaR);
  if (!full && !smooth) return out;
  NumericVector zin(N), zdg(N), vv(N), SS(N), ld(N);
  for (int k = 0; k < N; ++k) { zin[k] = vk[k] / std::sqrt(S[k] * sig2); zdg[k] = vk[k] / std::sqrt(Sung[k] * sig2); vv[k] = vk[k]; SS[k] = S[k] * sig2;
    ld[k] = -0.5 * (std::log(2.0 * M_PI * S[k] * sig2) + vk[k] * vk[k] / (S[k] * sig2)); }
  out["std_innov"] = zin; out["innovations_diagnostic"] = zdg; out["v"] = vv; out["S"] = SS; out["logdens"] = ld; out["gate_multiplier"] = gm_;
  Eigen::MatrixXd Ainv = Eigen::MatrixXd::Zero(m, m); bool ok_inv = A.allFinite() && std::fabs(A.determinant()) > 1e-300; if (ok_inv) Ainv = A.inverse();
  auto emit = [&](const std::vector<std::vector<V5> >& xs, const std::vector<M5>& Ps, const char* tag) {
    NumericVector xp(N), xsy(N), xo(N), xpv(N), xsv(N), xov(N), inet(N), xvv(N);
    for (int k = 0; k < N; ++k) {
      V5 c = xs[0][k]; for (int i = 0; i < m; ++i) c -= beta(i) * xs[i + 1][k];
      xp[k] = c(0); xsy[k] = c(1); xo[k] = c(2); inet[k] = c(4);
      double ex[3] = {0, 0, 0};
      if (ok_inv && sig2_fixed <= 0) for (int d = 0; d < 3; ++d) { Eigen::VectorXd xr(m); for (int i = 0; i < m; ++i) xr(i) = xs[i + 1][k](d); ex[d] = xr.dot(Ainv * xr); }
      xpv[k] = (Ps[k](0, 0) + ex[0]) * sig2; xsv[k] = (Ps[k](1, 1) + ex[1]) * sig2; xov[k] = (Ps[k](2, 2) + ex[2]) * sig2;
      double exv = 0; if (ok_inv && sig2_fixed <= 0) { Eigen::VectorXd xr(m); for (int i = 0; i < m; ++i) xr(i) = xs[i + 1][k](0) + xs[i + 1][k](2); exv = xr.dot(Ainv * xr); }
      xvv[k] = (Ps[k](0, 0) + Ps[k](2, 2) + 2.0 * Ps[k](0, 2) + exv) * sig2;                        // variance of the vagal total x_p + z1
    }
    std::string t(tag);
    out["X_p" + t] = xp; out["X_s" + t] = xsy; out["X_o" + t] = xo; out["X_p_var" + t] = xpv; out["X_s_var" + t] = xsv; out["X_o_var" + t] = xov; out["X_v_var" + t] = xvv; out["I_net" + t] = inet;
  };
  if (full) emit(xfilt, Pfilt, "");
  if (smooth) {                                                             // Rauch-Tung-Striebel on every stream (same gains)
    std::vector<std::vector<V5> > xs = xfilt; std::vector<M5> Ps = Pfilt;
    for (int k = N - 2; k >= 0; --k) {
      M5 G = Pfilt[k] * Fk[k + 1].transpose() * pinv5(Ppred[k + 1]);
      for (int j = 0; j < Mst; ++j) xs[j][k] = xfilt[j][k] + G * (xs[j][k + 1] - xpred[j][k + 1]);
      Ps[k] = Pfilt[k] + G * (Ps[k + 1] - Ppred[k + 1]) * G.transpose();
    }
    emit(xs, Ps, "_smooth");
  }
  return out;
}
