// EXPERIMENTAL (round 2, file reg3_filter.cpp). Regime-switching with regime-specific drive rates: regime j has its own (kappa_p,j, kappa_s,j) as well as targets mu_j; every regime keeps its own
// predicted covariance (full GPB1: per-regime Kalman branch, mixture likelihood, collapse after every beat). Optional shared damped oscillator. Derived from reg2_filter.cpp; the original header follows.
// Regime-switching targets (see ../reg_filter.cpp) with the optional damped oscillator of the M3O model (see ../../m3plus/ext_filter.cpp) and an approximate RTS smoother.
// State [x_p, x_s, z1, z2, I_net]; the regime j shifts the predicted mean by b_j = [(1-e_p) mu_p, (1-e_s) mu_s, 0, 0, -(dt-phi_p) mu_p + (dt-phi_s) mu_s]; covariances are regime independent.
// Forward pass: Gaussian mixture filter with collapse after every beat (GPB1). Smoother: Rauch-Tung-Striebel on the collapsed sequence, with the prior mixture mean and covariance as the
// one-step prediction (an approximation: it ignores that regimes are re-weighted by later data in the states, but the regime posterior W is returned filtered).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <complex>
#include <Eigen/SVD>
using namespace Rcpp;
typedef Eigen::Matrix<double, 5, 5> M5; typedef Eigen::Matrix<double, 5, 1> V5;
static M5 pinv5(const M5& X) {
  Eigen::JacobiSVD<M5> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV); V5 sv = svd.singularValues(); double tol = 1e-10 * std::max(sv(0), 1e-300); M5 Si = M5::Zero();
  for (int i = 0; i < 5; ++i) if (sv(i) > tol) Si(i, i) = 1.0 / sv(i);
  return svd.matrixV() * Si * svd.matrixU().transpose();
}
struct OU { double phi, q11, q13, q33; };
static inline OU ou(double k, double dt, double V) {
  double x = k * dt, s2 = 2.0 * k * V; OU b;
  if (x < 1e-3) { b.phi = dt - (k * dt * dt) / 2.0 + (k * k * dt * dt * dt) / 6.0; b.q11 = s2 * (dt - k * dt * dt + (2.0 * k * k * std::pow(dt, 3)) / 3.0);
    b.q13 = s2 * ((dt * dt) / 2.0 - (k * std::pow(dt, 3)) / 2.0 + (7.0 * k * k * std::pow(dt, 4)) / 24.0); b.q33 = s2 * (std::pow(dt, 3) / 3.0 - (k * std::pow(dt, 4)) / 4.0 + (7.0 * k * k * std::pow(dt, 5)) / 60.0); }
  else { double e1 = std::exp(-x), e2 = std::exp(-2.0 * x); b.phi = (1.0 - e1) / k; b.q11 = s2 * (1.0 - e2) / (2.0 * k); b.q13 = s2 * (1.0 - 2.0 * e1 + e2) / (2.0 * k * k); b.q33 = s2 * (2.0 * x - 3.0 + 4.0 * e1 - e2) / (2.0 * k * k * k); }
  return b;
}
struct Osc { double ea, c, s, gre, gim, A, Bre, Bim, C; };
static inline Osc osc_block(double alpha, double omega, double h, double Vo) {
  typedef std::complex<double> cd; cd lam(-alpha, omega), e = std::exp(lam * h), g = (e - 1.0) / lam; double base = (1.0 - std::exp(-2.0 * alpha * h)) / (2.0 * alpha), s2 = 4.0 * alpha * Vo; cd B = s2 / std::conj(lam) * (base - g);
  Osc o; o.ea = std::exp(-alpha * h); o.c = std::cos(omega * h); o.s = std::sin(omega * h); o.gre = g.real(); o.gim = g.imag(); o.A = s2 * base; o.Bre = B.real(); o.Bim = B.imag(); o.C = s2 / std::norm(lam) * (base - 2.0 * g.real() + h);
  return o;
}

// [[Rcpp::export]]
List reg3_filter(NumericVector dy_, NumericVector gm_, NumericVector kp_, NumericVector ks_, double lp, double lR, double sig2, double nu0, NumericMatrix mu, NumericMatrix Pi, NumericVector w0,
                 double osc_alpha, double osc_omega, double osc_lo, double lik_gate_max, bool full) {
  const int N = dy_.size(), K = mu.nrow(); const bool use_osc = osc_lo > 0 && osc_alpha > 0;
  V5 X = V5::Zero(); M5 P = M5::Zero(); P(0, 0) = lp * sig2; P(1, 1) = sig2; if (use_osc) { P(2, 2) = osc_lo * sig2; P(3, 3) = osc_lo * sig2; }
  std::vector<double> w(K); for (int j = 0; j < K; ++j) w[j] = w0[j];
  double LL = 0; NumericVector ld(N), zdg(N), xp(N), xs(N), xo(N), xvv(N), xpv(N), xsv(N); NumericMatrix W(N, K), Aprior(N, K);
  M5 I5 = M5::Identity(); std::vector<double> a(K), lj(K), wn(K), Sj(K); std::vector<V5> mj(K), Kj(K), mm(K); std::vector<M5> Ppj(K);
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], R = std::min(lR * gm_[k], 1e12) * dt * sig2, Rung = lR * dt * sig2, y = 1.0 - nu0 * dt;
    Osc o; if (use_osc) o = osc_block(osc_alpha, osc_omega, dt, osc_lo * sig2);
    double ymean = 0, mx = -1e300, Sung_mix = 0;
    for (int j = 0; j < K; ++j) {
      double s = 0; for (int i = 0; i < K; ++i) s += w[i] * Pi(i, j); a[j] = s; Aprior(k, j) = s;
      double kp = kp_[j], ks = ks_[j]; OU bp = ou(kp, dt, lp * sig2), bs = ou(ks, dt, sig2); double ep = std::exp(-kp * dt), es = std::exp(-ks * dt);
      M5 F = M5::Zero(), Q = M5::Zero(); F(0, 0) = ep; F(1, 1) = es; F(4, 0) = -bp.phi; F(4, 1) = bs.phi;
      Q(0, 0) = bp.q11; Q(0, 4) = Q(4, 0) = -bp.q13; Q(1, 1) = bs.q11; Q(1, 4) = Q(4, 1) = bs.q13; Q(4, 4) = bp.q33 + bs.q33;
      if (use_osc) { F(2, 2) = o.ea * o.c; F(2, 3) = -o.ea * o.s; F(3, 2) = o.ea * o.s; F(3, 3) = o.ea * o.c; F(4, 2) = -o.gre; F(4, 3) = o.gim;
        Q(2, 2) = o.A / 2; Q(3, 3) = o.A / 2; Q(2, 4) = Q(4, 2) = -o.Bre / 2; Q(3, 4) = Q(4, 3) = -o.Bim / 2; Q(4, 4) += o.C / 2; }
      M5 Pp = F * P * F.transpose() + Q; Ppj[j] = Pp; Sj[j] = std::max(Pp(4, 4) + R, 1e-12); Kj[j] = Pp.col(4) / Sj[j];
      V5 b = V5::Zero(); b(0) = (1.0 - ep) * mu(j, 0); b(1) = (1.0 - es) * mu(j, 1); b(4) = -(dt - bp.phi) * mu(j, 0) + (dt - bs.phi) * mu(j, 1); mj[j] = F * X + b;
      double v = y - mj[j](4); lj[j] = std::log(std::max(a[j], 1e-300)) - 0.5 * v * v / Sj[j] - 0.5 * std::log(2.0 * M_PI * Sj[j]); if (lj[j] > mx) mx = lj[j];
      ymean += a[j] * mj[j](4); Sung_mix += a[j] * std::max(Pp(4, 4) + Rung, 1e-12);
    }
    double sum = 0; for (int j = 0; j < K; ++j) sum += std::exp(lj[j] - mx); double logp = mx + std::log(sum); ld[k] = logp; if (gm_[k] < lik_gate_max) LL += logp;
    for (int j = 0; j < K; ++j) wn[j] = std::exp(lj[j] - logp);
    double spread = 0; for (int j = 0; j < K; ++j) spread += a[j] * (mj[j](4) - ymean) * (mj[j](4) - ymean); zdg[k] = (y - ymean) / std::sqrt(Sung_mix + spread);
    V5 Xn = V5::Zero(); for (int j = 0; j < K; ++j) { mm[j] = mj[j] + Kj[j] * (y - mj[j](4)); Xn += wn[j] * mm[j]; }
    M5 Pn = M5::Zero(); for (int j = 0; j < K; ++j) { M5 IKH = I5; IKH.col(4) -= Kj[j]; M5 Pj = IKH * Ppj[j] * IKH.transpose() + Kj[j] * Kj[j].transpose() * R; V5 d = mm[j] - Xn; Pn += wn[j] * (Pj + d * d.transpose()); }
    X = Xn; P = 0.5 * (Pn + Pn.transpose()); for (int j = 0; j < K; ++j) { w[j] = wn[j]; W(k, j) = wn[j]; }
    xp[k] = X(0); xs[k] = X(1); xo[k] = X(2); xpv[k] = P(0, 0); xsv[k] = P(1, 1); xvv[k] = P(0, 0) + P(2, 2) + 2.0 * P(0, 2);
  }
  List out = List::create(_["LL"] = LL); if (!full) return out;
  out["logdens"] = ld; out["z"] = zdg; out["W"] = W; out["A"] = Aprior; out["X_p"] = xp; out["X_s"] = xs; out["X_o"] = xo; out["X_p_var"] = xpv; out["X_s_var"] = xsv; out["X_v_var"] = xvv;
  return out;
}
