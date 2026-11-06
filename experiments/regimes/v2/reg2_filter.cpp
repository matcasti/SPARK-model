// EXPERIMENTAL (round 2). Regime-switching targets (see ../reg_filter.cpp) with the optional damped oscillator of the M3O model (see ../../m3plus/ext_filter.cpp) and an approximate RTS smoother.
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
List reg2_filter(NumericVector dy_, NumericVector gm_, double kp, double ks, double lp, double lR, double sig2, double nu0, NumericMatrix mu, NumericMatrix Pi, NumericVector w0,
                 double osc_alpha, double osc_omega, double osc_lo, double lik_gate_max, bool full, bool smooth, double eps = 0.0, double kappa = 25.0) {
  const int N = dy_.size(), K = mu.nrow(); const bool use_osc = osc_lo > 0 && osc_alpha > 0;
  V5 X = V5::Zero(); M5 P = M5::Zero(); P(0, 0) = lp * sig2; P(1, 1) = sig2; if (use_osc) { P(2, 2) = osc_lo * sig2; P(3, 3) = osc_lo * sig2; }
  std::vector<double> w(K); for (int j = 0; j < K; ++j) w[j] = w0[j];
  double LL = 0; NumericVector ld(N), ldr(N), zdg(N), xp(N), xs(N), xo(N), xvv(N), xpv(N), xsv(N); NumericMatrix W(N, K), Aprior(N, K);
  std::vector<V5> Xf, Xpr; std::vector<M5> Pf, Ppr, Fk; if (smooth) { Xf.resize(N); Xpr.resize(N); Pf.resize(N); Ppr.resize(N); Fk.resize(N); }
  M5 I5 = M5::Identity(); std::vector<double> a(K), lj(K), wn(K); std::vector<V5> mj(K), bj(K);
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], R = std::min(lR * gm_[k], 1e12) * dt * sig2, Rung = lR * dt * sig2;
    OU bp = ou(kp, dt, lp * sig2), bs = ou(ks, dt, sig2); double ep = std::exp(-kp * dt), es = std::exp(-ks * dt);
    M5 F = M5::Zero(), Q = M5::Zero(); F(0, 0) = ep; F(1, 1) = es; F(4, 0) = -bp.phi; F(4, 1) = bs.phi;
    Q(0, 0) = bp.q11; Q(0, 4) = Q(4, 0) = -bp.q13; Q(1, 1) = bs.q11; Q(1, 4) = Q(4, 1) = bs.q13; Q(4, 4) = bp.q33 + bs.q33;
    if (use_osc) { Osc o = osc_block(osc_alpha, osc_omega, dt, osc_lo * sig2);
      F(2, 2) = o.ea * o.c; F(2, 3) = -o.ea * o.s; F(3, 2) = o.ea * o.s; F(3, 3) = o.ea * o.c; F(4, 2) = -o.gre; F(4, 3) = o.gim;
      Q(2, 2) = o.A / 2; Q(3, 3) = o.A / 2; Q(2, 4) = Q(4, 2) = -o.Bre / 2; Q(3, 4) = Q(4, 3) = -o.Bim / 2; Q(4, 4) += o.C / 2; }
    M5 Pp = F * P * F.transpose() + Q; double S = std::max(Pp(4, 4) + R, 1e-12), Sung = std::max(Pp(4, 4) + Rung, 1e-12); V5 Kg = Pp.col(4) / S; double y = 1.0 - nu0 * dt;
    V5 Xp = F * X; double ymean = 0, mx = -1e300; V5 bbar = V5::Zero();
    for (int j = 0; j < K; ++j) { double s = 0; for (int i = 0; i < K; ++i) s += w[i] * Pi(i, j); a[j] = s; Aprior(k, j) = s;
      V5 b = V5::Zero(); b(0) = (1.0 - ep) * mu(j, 0); b(1) = (1.0 - es) * mu(j, 1); b(4) = -(dt - bp.phi) * mu(j, 0) + (dt - bs.phi) * mu(j, 1); bj[j] = b; bbar += a[j] * b; mj[j] = Xp + b; ymean += a[j] * mj[j](4);
      double v = y - mj[j](4); lj[j] = std::log(std::max(a[j], 1e-300)) - 0.5 * v * v / S - 0.5 * std::log(2.0 * M_PI * S); if (lj[j] > mx) mx = lj[j]; }
    double sum = 0; for (int j = 0; j < K; ++j) sum += std::exp(lj[j] - mx); double logp = mx + std::log(sum); ld[k] = logp; if (gm_[k] < lik_gate_max) LL += logp;
    for (int j = 0; j < K; ++j) wn[j] = std::exp(lj[j] - logp);
    double spread_pred = 0; M5 Ptot = Pp; for (int j = 0; j < K; ++j) { spread_pred += a[j] * (mj[j](4) - ymean) * (mj[j](4) - ymean); V5 d = bj[j] - bbar; Ptot += a[j] * d * d.transpose(); }
    zdg[k] = (y - ymean) / std::sqrt(Sung + spread_pred);
    { double Vm = S + spread_pred, vv = y - ymean, lo = -0.5 * vv * vv / (kappa * Vm) - 0.5 * std::log(2.0 * M_PI * kappa * Vm); if (eps > 0) { double a1 = std::log1p(-eps) + logp, a2 = std::log(eps) + lo, m = std::max(a1, a2); ldr[k] = m + std::log(std::exp(a1 - m) + std::exp(a2 - m)); } else ldr[k] = logp; }
    if (smooth) { Xpr[k] = Xp + bbar; Ppr[k] = Ptot; Fk[k] = F; }
    V5 Xn = V5::Zero(); std::vector<V5> mm(K); for (int j = 0; j < K; ++j) { mm[j] = mj[j] + Kg * (y - mj[j](4)); Xn += wn[j] * mm[j]; }
    M5 IKH = I5; IKH.col(4) -= Kg; M5 Pn = IKH * Pp * IKH.transpose() + Kg * Kg.transpose() * R; for (int j = 0; j < K; ++j) { V5 d = mm[j] - Xn; Pn += wn[j] * d * d.transpose(); }
    X = Xn; P = 0.5 * (Pn + Pn.transpose()); for (int j = 0; j < K; ++j) { w[j] = wn[j]; W(k, j) = wn[j]; }
    if (smooth) { Xf[k] = X; Pf[k] = P; }
    xp[k] = X(0); xs[k] = X(1); xo[k] = X(2); xpv[k] = P(0, 0); xsv[k] = P(1, 1); xvv[k] = P(0, 0) + P(2, 2) + 2.0 * P(0, 2);
  }
  List out = List::create(_["LL"] = LL); if (!full && !smooth) return out;
  out["logdens"] = ld; out["logdens_rob"] = ldr; out["z"] = zdg; out["W"] = W; out["A"] = Aprior; out["X_p"] = xp; out["X_s"] = xs; out["X_o"] = xo; out["X_p_var"] = xpv; out["X_s_var"] = xsv; out["X_v_var"] = xvv;
  if (smooth) {
    std::vector<V5> Xs = Xf; std::vector<M5> Ps = Pf;
    for (int k = N - 2; k >= 0; --k) { M5 G = Pf[k] * Fk[k + 1].transpose() * pinv5(Ppr[k + 1]); Xs[k] = Xf[k] + G * (Xs[k + 1] - Xpr[k + 1]); Ps[k] = Pf[k] + G * (Ps[k + 1] - Ppr[k + 1]) * G.transpose(); }
    NumericVector sp(N), ss(N), so(N), spv(N), ssv(N), svv(N);
    for (int k = 0; k < N; ++k) { sp[k] = Xs[k](0); ss[k] = Xs[k](1); so[k] = Xs[k](2); spv[k] = Ps[k](0, 0); ssv[k] = Ps[k](1, 1); svv[k] = Ps[k](0, 0) + Ps[k](2, 2) + 2.0 * Ps[k](0, 2); }
    out["Xs_p"] = sp; out["Xs_s"] = ss; out["Xs_o"] = so; out["Xs_p_var"] = spv; out["Xs_s_var"] = ssv; out["Xs_v_var"] = svv;
  }
  return out;
}

// Smoother with soft regime inputs: given smoothed regime probabilities Ws (N x K; e.g. from Kim's backward recursion on the filtered W), run a linear Gaussian filter whose predicted mean is shifted
// by sum_j Ws[k, j] b_j (the prior variance is inflated by the spread of the b_j), then RTS. States only; no likelihood.
// [[Rcpp::export]]
List reg2_smooth_states(NumericVector dy_, NumericVector gm_, double kp, double ks, double lp, double lR, double sig2, double nu0, NumericMatrix mu, NumericMatrix Ws,
                        double osc_alpha, double osc_omega, double osc_lo) {
  const int N = dy_.size(), K = mu.nrow(); const bool use_osc = osc_lo > 0 && osc_alpha > 0;
  V5 X = V5::Zero(); M5 P = M5::Zero(); P(0, 0) = lp * sig2; P(1, 1) = sig2; if (use_osc) { P(2, 2) = osc_lo * sig2; P(3, 3) = osc_lo * sig2; }
  std::vector<V5> Xf(N), Xpr(N); std::vector<M5> Pf(N), Ppr(N), Fk(N); M5 I5 = M5::Identity();
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], R = std::min(lR * gm_[k], 1e12) * dt * sig2; OU bp = ou(kp, dt, lp * sig2), bs = ou(ks, dt, sig2); double ep = std::exp(-kp * dt), es = std::exp(-ks * dt);
    M5 F = M5::Zero(), Q = M5::Zero(); F(0, 0) = ep; F(1, 1) = es; F(4, 0) = -bp.phi; F(4, 1) = bs.phi;
    Q(0, 0) = bp.q11; Q(0, 4) = Q(4, 0) = -bp.q13; Q(1, 1) = bs.q11; Q(1, 4) = Q(4, 1) = bs.q13; Q(4, 4) = bp.q33 + bs.q33;
    if (use_osc) { Osc o = osc_block(osc_alpha, osc_omega, dt, osc_lo * sig2);
      F(2, 2) = o.ea * o.c; F(2, 3) = -o.ea * o.s; F(3, 2) = o.ea * o.s; F(3, 3) = o.ea * o.c; F(4, 2) = -o.gre; F(4, 3) = o.gim;
      Q(2, 2) = o.A / 2; Q(3, 3) = o.A / 2; Q(2, 4) = Q(4, 2) = -o.Bre / 2; Q(3, 4) = Q(4, 3) = -o.Bim / 2; Q(4, 4) += o.C / 2; }
    std::vector<V5> bj(K); V5 bbar = V5::Zero();
    for (int j = 0; j < K; ++j) { V5 b = V5::Zero(); b(0) = (1.0 - ep) * mu(j, 0); b(1) = (1.0 - es) * mu(j, 1); b(4) = -(dt - bp.phi) * mu(j, 0) + (dt - bs.phi) * mu(j, 1); bj[j] = b; bbar += Ws(k, j) * b; }
    M5 Pp = F * P * F.transpose() + Q; for (int j = 0; j < K; ++j) { V5 d = bj[j] - bbar; Pp += Ws(k, j) * d * d.transpose(); }
    double S = std::max(Pp(4, 4) + R, 1e-12); V5 Kg = Pp.col(4) / S; double y = 1.0 - nu0 * dt; V5 xp = F * X + bbar; X = xp + Kg * (y - xp(4));
    M5 IKH = I5; IKH.col(4) -= Kg; P = IKH * Pp * IKH.transpose() + Kg * Kg.transpose() * R; P = 0.5 * (P + P.transpose()); Xf[k] = X; Pf[k] = P; Xpr[k] = xp; Ppr[k] = Pp; Fk[k] = F;
  }
  std::vector<V5> Xs = Xf; std::vector<M5> Ps = Pf;
  for (int k = N - 2; k >= 0; --k) { M5 G = Pf[k] * Fk[k + 1].transpose() * pinv5(Ppr[k + 1]); Xs[k] = Xf[k] + G * (Xs[k + 1] - Xpr[k + 1]); Ps[k] = Pf[k] + G * (Ps[k + 1] - Ppr[k + 1]) * G.transpose(); }
  NumericVector sp(N), ss(N), so(N), spv(N), ssv(N), svv(N);
  for (int k = 0; k < N; ++k) { sp[k] = Xs[k](0); ss[k] = Xs[k](1); so[k] = Xs[k](2); spv[k] = Ps[k](0, 0); ssv[k] = Ps[k](1, 1); svv[k] = Ps[k](0, 0) + Ps[k](2, 2) + 2.0 * Ps[k](0, 2); }
  return List::create(_["X_p"] = sp, _["X_s"] = ss, _["X_o"] = so, _["X_p_var"] = spv, _["X_s_var"] = ssv, _["X_v_var"] = svv);
}
