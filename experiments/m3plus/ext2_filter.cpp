// EXPERIMENTAL filter with stochastic SET-POINTS: each drive relaxes (rate kappa) towards a target m_b that moves as a random walk (dm = s_m dW), so a change of the
// target is a "stochastic step" for the fast vagal drive and a slow "ramp" for the sympathetic one: dx = kappa (m - x) dt + s_x dW,  I = integral of x.
// State S = [x_p, x_s, z1, z2, m_p, m_s, I_net]; optional damped oscillator (z1, z2) as in ext_filter.cpp. Process covariances by 10-point Gauss-Legendre quadrature.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <complex>
#include <Eigen/SVD>
using namespace Rcpp;
typedef Eigen::Matrix<double, 7, 7> M7; typedef Eigen::Matrix<double, 7, 1> V7;
static M7 pinv7(const M7& X) { Eigen::JacobiSVD<M7> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV); V7 sv = svd.singularValues(); double tol = 1e-10 * std::max(sv(0), 1e-300); M7 Si = M7::Zero();
  for (int i = 0; i < 7; ++i) if (sv(i) > tol) Si(i, i) = 1.0 / sv(i); return svd.matrixV() * Si * svd.matrixU().transpose(); }
static const double GLX[10] = {0.013046735741414, 0.067468316655507, 0.160295215850488, 0.283302302935376, 0.425562830509184, 0.574437169490816, 0.716697697064624, 0.839704784149512, 0.932531683344493, 0.986953264258586};
static const double GLW[10] = {0.033335672154344, 0.074725674575290, 0.109543181257991, 0.134633359654998, 0.147762112357376, 0.147762112357376, 0.134633359654998, 0.109543181257991, 0.074725674575290, 0.033335672154344};
// Block of one drive: noise intensities sx2 (OU), sm2 (set-point random walk). Output: transition entries and process covariance of (x, m, I).
struct SP { double e, phi, h; double vxx, vxi, vii, vxm, vim, vmm; };
static inline SP sp_block(double kappa, double h, double sx2, double sm2) {
  SP o; o.h = h; o.e = std::exp(-kappa * h); o.phi = (1.0 - o.e) / kappa; o.vxx = o.vxi = o.vii = o.vxm = o.vim = 0; o.vmm = sm2 * h;
  for (int i = 0; i < 10; ++i) { double t = GLX[i] * h, w = GLW[i] * h, e = std::exp(-kappa * t), a = 1.0 - e, ph = a / kappa, b = t - ph;
    o.vxx += w * (sx2 * e * e + sm2 * a * a); o.vxi += w * (sx2 * e * ph + sm2 * a * b); o.vii += w * (sx2 * ph * ph + sm2 * b * b); o.vxm += w * sm2 * a; o.vim += w * sm2 * b; }
  return o;
}
//' Set-point block for tests: rows/cols x, m, I.
// [[Rcpp::export]]
List sp_moments(double kappa, double h, double sx2, double sm2) { SP o = sp_block(kappa, h, sx2, sm2); NumericMatrix F(3, 3), Q(3, 3);
  F(0, 0) = o.e; F(0, 1) = 1 - o.e; F(1, 1) = 1; F(2, 0) = o.phi; F(2, 1) = h - o.phi;
  Q(0, 0) = o.vxx; Q(0, 1) = Q(1, 0) = o.vxm; Q(0, 2) = Q(2, 0) = o.vxi; Q(1, 1) = o.vmm; Q(1, 2) = Q(2, 1) = o.vim; Q(2, 2) = o.vii; return List::create(_["F"] = F, _["Q"] = Q); }
struct Osc { double ea, c, s, gre, gim, A, Bre, Bim, C; };
static inline Osc osc_block(double alpha, double omega, double h, double Vo) {
  typedef std::complex<double> cd; cd lam(-alpha, omega), e = std::exp(lam * h), g = (e - 1.0) / lam; double base = (1.0 - std::exp(-2.0 * alpha * h)) / (2.0 * alpha), s2 = 4.0 * alpha * Vo; cd B = s2 / std::conj(lam) * (base - g);
  Osc o; o.ea = std::exp(-alpha * h); o.c = std::cos(omega * h); o.s = std::sin(omega * h); o.gre = g.real(); o.gim = g.imag(); o.A = s2 * base; o.Bre = B.real(); o.Bim = B.imag(); o.C = s2 / std::norm(lam) * (base - 2.0 * g.real() + h); return o; }

// [[Rcpp::export]]
List ext2_filter(NumericVector dy_, NumericVector gm_, NumericVector qs_, NumericMatrix Xreg, double kp, double ks, double lp, double lR,
                 double osc_alpha, double osc_omega, double osc_lo, double q_p, double q_s, double lik_gate_max, bool full, bool smooth, NumericVector beta_fixed, double sig2_fixed) {
  const int N = dy_.size(), m = Xreg.ncol(), Mst = m + 1; const bool use_osc = osc_lo > 0 && osc_alpha > 0;
  std::vector<V7> Xf(Mst, V7::Zero()); M7 P = M7::Zero(); P(0, 0) = lp; P(1, 1) = 1.0; if (use_osc) { P(2, 2) = osc_lo; P(3, 3) = osc_lo; }
  std::vector<double> S(N), Sung(N); std::vector<std::vector<double> > v(Mst, std::vector<double>(N)); bool keep = full || smooth;
  std::vector<std::vector<V7> > xfilt, xpred; std::vector<M7> Pfilt, Ppred, Fk; if (keep) { xfilt.assign(Mst, std::vector<V7>(N)); xpred.assign(Mst, std::vector<V7>(N)); Pfilt.resize(N); Ppred.resize(N); Fk.resize(N); }
  M7 I7 = M7::Identity();
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], qs = qs_[k], R = std::min(lR * gm_[k], 1e12) * dt;
    SP bp = sp_block(kp, dt, 2.0 * kp * lp * qs, q_p * qs), bs = sp_block(ks, dt, 2.0 * ks * 1.0 * qs, q_s * qs);
    M7 F = M7::Zero(), Q = M7::Zero();
    F(0, 0) = bp.e; F(0, 4) = 1.0 - bp.e; F(4, 4) = 1.0; F(1, 1) = bs.e; F(1, 5) = 1.0 - bs.e; F(5, 5) = 1.0;
    F(6, 0) = -bp.phi; F(6, 4) = -(dt - bp.phi); F(6, 1) = bs.phi; F(6, 5) = dt - bs.phi;
    Q(0, 0) = bp.vxx; Q(0, 4) = Q(4, 0) = bp.vxm; Q(4, 4) = bp.vmm; Q(0, 6) = Q(6, 0) = -bp.vxi; Q(4, 6) = Q(6, 4) = -bp.vim;
    Q(1, 1) = bs.vxx; Q(1, 5) = Q(5, 1) = bs.vxm; Q(5, 5) = bs.vmm; Q(1, 6) = Q(6, 1) = bs.vxi; Q(5, 6) = Q(6, 5) = bs.vim; Q(6, 6) = bp.vii + bs.vii;
    if (use_osc) { Osc o = osc_block(osc_alpha, osc_omega, dt, osc_lo); F(2, 2) = o.ea * o.c; F(2, 3) = -o.ea * o.s; F(3, 2) = o.ea * o.s; F(3, 3) = o.ea * o.c; F(6, 2) = -o.gre; F(6, 3) = o.gim;
      Q(2, 2) = o.A / 2 * qs; Q(3, 3) = o.A / 2 * qs; Q(2, 6) = Q(6, 2) = -o.Bre / 2 * qs; Q(3, 6) = Q(6, 3) = -o.Bim / 2 * qs; Q(6, 6) += o.C / 2 * qs; }
    M7 Pp = F * P * F.transpose() + Q; double Sk = std::max(Pp(6, 6) + R, 1e-12); S[k] = Sk; Sung[k] = std::max(Pp(6, 6) + lR * dt, 1e-12); V7 K = Pp.col(6) / Sk;
    for (int j = 0; j < Mst; ++j) { V7 xp = F * Xf[j]; double target = (j == 0) ? 1.0 : Xreg(k, j - 1); double vj = target - xp(6); v[j][k] = vj; if (keep) xpred[j][k] = xp; Xf[j] = xp + K * vj; if (keep) xfilt[j][k] = Xf[j]; }
    M7 IKH = I7; IKH.col(6) -= K; P = IKH * Pp * IKH.transpose() + K * K.transpose() * R; P = 0.5 * (P + P.transpose()); if (keep) { Pfilt[k] = P; Ppred[k] = Pp; Fk[k] = F; }
  }
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); Eigen::VectorXd b = Eigen::VectorXd::Zero(m);
  for (int k = 0; k < N; ++k) for (int i = 0; i < m; ++i) { b(i) += v[i + 1][k] * v[0][k] / S[k]; for (int j = 0; j < m; ++j) A(i, j) += v[i + 1][k] * v[j + 1][k] / S[k]; }
  Eigen::VectorXd beta(m);
  if (beta_fixed.size() == m) { for (int i = 0; i < m; ++i) beta(i) = beta_fixed[i]; }
  else { bool fin = A.allFinite() && b.allFinite(); if (fin) beta = A.colPivHouseholderQr().solve(b); if (!fin || !beta.allFinite() || !(beta(0) > 0)) return List::create(_["LL_star"] = R_NegInf, _["nu0"] = NA_REAL, _["sig2"] = NA_REAL); }
  std::vector<double> vk(N); int n_eff = 0; for (int k = 0; k < N; ++k) if (gm_[k] < lik_gate_max) ++n_eff; bool all = n_eff < 2; if (all) n_eff = N; double s2 = 0, sumlogS = 0;
  for (int k = 0; k < N; ++k) { double r = v[0][k]; for (int i = 0; i < m; ++i) r -= beta(i) * v[i + 1][k]; vk[k] = r; if (all || gm_[k] < lik_gate_max) { s2 += r * r / S[k]; sumlogS += std::log(S[k]); } }
  double sig2 = (sig2_fixed > 0) ? sig2_fixed : std::max(s2 / n_eff, 1e-12); double LL = -0.5 * sumlogS - (n_eff / 2.0) * std::log(sig2);
  NumericVector betaR(m); for (int i = 0; i < m; ++i) betaR[i] = beta(i);
  List out = List::create(_["LL_star"] = LL, _["nu0"] = beta(0), _["sig2"] = sig2, _["beta"] = betaR); if (!full && !smooth) return out;
  NumericVector zin(N), zdg(N), vv(N), SS(N), ld(N);
  for (int k = 0; k < N; ++k) { zin[k] = vk[k] / std::sqrt(S[k] * sig2); zdg[k] = vk[k] / std::sqrt(Sung[k] * sig2); vv[k] = vk[k]; SS[k] = S[k] * sig2; ld[k] = -0.5 * (std::log(2.0 * M_PI * S[k] * sig2) + vk[k] * vk[k] / (S[k] * sig2)); }
  out["std_innov"] = zin; out["innovations_diagnostic"] = zdg; out["v"] = vv; out["S"] = SS; out["logdens"] = ld; out["gate_multiplier"] = gm_;
  auto emit = [&](const std::vector<std::vector<V7> >& xs, const std::vector<M7>& Ps, const char* tag) {
    NumericVector xp(N), xsy(N), xo(N), mp(N), ms(N), xpv(N), xsv(N), xov(N), xvv(N), mpv(N), msv(N);
    for (int k = 0; k < N; ++k) { V7 c = xs[0][k]; for (int i = 0; i < m; ++i) c -= beta(i) * xs[i + 1][k]; xp[k] = c(0); xsy[k] = c(1); xo[k] = c(2); mp[k] = c(4); ms[k] = c(5);
      xpv[k] = Ps[k](0, 0) * sig2; xsv[k] = Ps[k](1, 1) * sig2; xov[k] = Ps[k](2, 2) * sig2; xvv[k] = (Ps[k](0, 0) + Ps[k](2, 2) + 2.0 * Ps[k](0, 2)) * sig2; mpv[k] = Ps[k](4, 4) * sig2; msv[k] = Ps[k](5, 5) * sig2; }
    std::string t(tag); out["X_p" + t] = xp; out["X_s" + t] = xsy; out["X_o" + t] = xo; out["M_p" + t] = mp; out["M_s" + t] = ms; out["X_p_var" + t] = xpv; out["X_s_var" + t] = xsv; out["X_o_var" + t] = xov; out["X_v_var" + t] = xvv; out["M_p_var" + t] = mpv; out["M_s_var" + t] = msv; };
  if (full) emit(xfilt, Pfilt, "");
  if (smooth) { std::vector<std::vector<V7> > xs = xfilt; std::vector<M7> Ps = Pfilt;
    for (int k = N - 2; k >= 0; --k) { M7 G = Pfilt[k] * Fk[k + 1].transpose() * pinv7(Ppred[k + 1]); for (int j = 0; j < Mst; ++j) xs[j][k] = xfilt[j][k] + G * (xs[j][k + 1] - xpred[j][k + 1]); Ps[k] = Pfilt[k] + G * (Ps[k + 1] - Ppred[k + 1]) * G.transpose(); }
    emit(xs, Ps, "_smooth"); }
  return out;
}
