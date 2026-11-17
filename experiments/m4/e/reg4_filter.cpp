// EXPERIMENTAL (M4e, E2). Regime-switching set-points with their own response kinetics. State [u_p, u_s, z1, z2, I_net, m_p, m_s]:
//   u_b  OU deviation of drive b (rate kappa_b, as in M3); z1, z2 the optional damped respiratory oscillator (vagal pathway);
//   m_b  deterministic set-point response: m_b relaxes towards mu_b(s) at rate rate(s, b), s = the regime of the current beat interval (so rate(1, b) = "on", rate(0, b) = "off").
//   Vagal-like state = u_p + m_p (+ z1), sympathetic-like = u_s + m_s; I_net integrates -(vagal) + sympathetic; y_k = 1 - nu0 dt_k = I_net + eta_k.
// With rate(s, b) = kappa_b for every s this is the model of ../../regimes/v2/reg2_filter.cpp. Forward pass: GPB1 with regime-specific transition matrices (collapse after every beat).
// Smoother: RTS on the collapsed sequence, with the prior-weighted transition matrix and the predicted mixture moments (an approximation).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <complex>
#include <Eigen/SVD>
using namespace Rcpp;
typedef Eigen::Matrix<double, 7, 7> M7; typedef Eigen::Matrix<double, 7, 1> V7;
static M7 pinv7(const M7& X) {
  Eigen::JacobiSVD<M7> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV); V7 sv = svd.singularValues(); double tol = 1e-10 * std::max(sv(0), 1e-300); M7 Si = M7::Zero();
  for (int i = 0; i < 7; ++i) if (sv(i) > tol) Si(i, i) = 1.0 / sv(i);
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
static inline double phi_of(double r, double dt) { double x = r * dt; return x < 1e-6 ? dt - r * dt * dt / 2.0 : (1.0 - std::exp(-x)) / r; }
struct Osc { double ea, c, s, gre, gim, A, Bre, Bim, C; };
static inline Osc osc_block(double alpha, double omega, double h, double Vo) {
  typedef std::complex<double> cd; cd lam(-alpha, omega), e = std::exp(lam * h), g = (e - 1.0) / lam; double base = (1.0 - std::exp(-2.0 * alpha * h)) / (2.0 * alpha), s2 = 4.0 * alpha * Vo; cd B = s2 / std::conj(lam) * (base - g);
  Osc o; o.ea = std::exp(-alpha * h); o.c = std::cos(omega * h); o.s = std::sin(omega * h); o.gre = g.real(); o.gim = g.imag(); o.A = s2 * base; o.Bre = B.real(); o.Bim = B.imag(); o.C = s2 / std::norm(lam) * (base - 2.0 * g.real() + h);
  return o;
}
// [[Rcpp::export]]
List reg4_filter(NumericVector dy_, NumericVector gm_, double kp, double ks, double lp, double lR, double sig2, double nu0, NumericMatrix mu, NumericMatrix rate, NumericMatrix Pi, NumericVector w0,
                 double osc_alpha, double osc_omega, double osc_lo, double lik_gate_max, bool full, bool smooth) {
  const int N = dy_.size(), K = mu.nrow(); const bool use_osc = osc_lo > 0 && osc_alpha > 0;
  V7 X = V7::Zero(); M7 P = M7::Zero(); P(0, 0) = lp * sig2; P(1, 1) = sig2; if (use_osc) { P(2, 2) = osc_lo * sig2; P(3, 3) = osc_lo * sig2; }
  std::vector<double> w(K); for (int j = 0; j < K; ++j) w[j] = w0[j];
  double LL = 0; NumericVector ld(N), zdg(N), xp(N), xs(N), xo(N), xvv(N), xpv(N), xsv(N), mp(N), ms(N); NumericMatrix W(N, K);
  std::vector<V7> Xf, Xpr; std::vector<M7> Pf, Ppr, Fk; if (smooth) { Xf.resize(N); Xpr.resize(N); Pf.resize(N); Ppr.resize(N); Fk.resize(N); }
  M7 I7 = M7::Identity(); std::vector<double> a(K), lj(K), wn(K), Sj(K), Suj(K); std::vector<V7> mj(K), mm(K); std::vector<M7> Pj(K), Pnj(K), Fj(K);
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], R = std::min(lR * gm_[k], 1e12) * dt * sig2, Rung = lR * dt * sig2; double y = 1.0 - nu0 * dt;
    OU bp = ou(kp, dt, lp * sig2), bs = ou(ks, dt, sig2); double ep = std::exp(-kp * dt), es = std::exp(-ks * dt);
    M7 F0 = M7::Zero(), Q = M7::Zero(); F0(0, 0) = ep; F0(1, 1) = es; F0(4, 0) = -bp.phi; F0(4, 1) = bs.phi;
    Q(0, 0) = bp.q11; Q(0, 4) = Q(4, 0) = -bp.q13; Q(1, 1) = bs.q11; Q(1, 4) = Q(4, 1) = bs.q13; Q(4, 4) = bp.q33 + bs.q33;
    if (use_osc) { Osc o = osc_block(osc_alpha, osc_omega, dt, osc_lo * sig2);
      F0(2, 2) = o.ea * o.c; F0(2, 3) = -o.ea * o.s; F0(3, 2) = o.ea * o.s; F0(3, 3) = o.ea * o.c; F0(4, 2) = -o.gre; F0(4, 3) = o.gim;
      Q(2, 2) = o.A / 2; Q(3, 3) = o.A / 2; Q(2, 4) = Q(4, 2) = -o.Bre / 2; Q(3, 4) = Q(4, 3) = -o.Bim / 2; Q(4, 4) += o.C / 2; }
    double mx = -1e300, ymean = 0;
    for (int j = 0; j < K; ++j) { double s = 0; for (int i = 0; i < K; ++i) s += w[i] * Pi(i, j); a[j] = s;
      double rp = rate(j, 0), rs = rate(j, 1), erp = std::exp(-rp * dt), ers = std::exp(-rs * dt), php = phi_of(rp, dt), phs = phi_of(rs, dt);
      M7 F = F0; F(5, 5) = erp; F(6, 6) = ers; F(4, 5) = -php; F(4, 6) = phs; Fj[j] = F;
      V7 b = V7::Zero(); b(5) = (1.0 - erp) * mu(j, 0); b(6) = (1.0 - ers) * mu(j, 1); b(4) = -(dt - php) * mu(j, 0) + (dt - phs) * mu(j, 1);
      mj[j] = F * X + b; Pj[j] = F * P * F.transpose() + Q; Sj[j] = std::max(Pj[j](4, 4) + R, 1e-12); Suj[j] = std::max(Pj[j](4, 4) + Rung, 1e-12); ymean += a[j] * mj[j](4);
      double v = y - mj[j](4); lj[j] = std::log(std::max(a[j], 1e-300)) - 0.5 * v * v / Sj[j] - 0.5 * std::log(2.0 * M_PI * Sj[j]); if (lj[j] > mx) mx = lj[j]; }
    double sum = 0; for (int j = 0; j < K; ++j) sum += std::exp(lj[j] - mx); double logp = mx + std::log(sum); ld[k] = logp; if (gm_[k] < lik_gate_max) LL += logp;
    double vpred = 0; for (int j = 0; j < K; ++j) { wn[j] = std::exp(lj[j] - logp); vpred += a[j] * (Suj[j] + (mj[j](4) - ymean) * (mj[j](4) - ymean)); }
    zdg[k] = (y - ymean) / std::sqrt(vpred);
    if (smooth) { V7 xb = V7::Zero(); M7 Fb = M7::Zero(); for (int j = 0; j < K; ++j) { xb += a[j] * mj[j]; Fb += a[j] * Fj[j]; }
      M7 Pb = M7::Zero(); for (int j = 0; j < K; ++j) { V7 d = mj[j] - xb; Pb += a[j] * (Pj[j] + d * d.transpose()); } Xpr[k] = xb; Ppr[k] = Pb; Fk[k] = Fb; }
    V7 Xn = V7::Zero();
    for (int j = 0; j < K; ++j) { V7 Kg = Pj[j].col(4) / Sj[j]; mm[j] = mj[j] + Kg * (y - mj[j](4)); M7 IKH = I7; IKH.col(4) -= Kg; Pnj[j] = IKH * Pj[j] * IKH.transpose() + Kg * Kg.transpose() * R; Xn += wn[j] * mm[j]; }
    M7 Pn = M7::Zero(); for (int j = 0; j < K; ++j) { V7 d = mm[j] - Xn; Pn += wn[j] * (Pnj[j] + d * d.transpose()); }
    X = Xn; P = 0.5 * (Pn + Pn.transpose()); for (int j = 0; j < K; ++j) { w[j] = wn[j]; W(k, j) = wn[j]; }
    if (smooth) { Xf[k] = X; Pf[k] = P; }
    xp[k] = X(0) + X(5); xs[k] = X(1) + X(6); xo[k] = X(2); mp[k] = X(5); ms[k] = X(6);
    xpv[k] = P(0, 0) + P(5, 5) + 2.0 * P(0, 5); xsv[k] = P(1, 1) + P(6, 6) + 2.0 * P(1, 6); xvv[k] = xpv[k] + P(2, 2) + 2.0 * (P(0, 2) + P(2, 5));
  }
  List out = List::create(_["LL"] = LL); if (!full && !smooth) return out;
  out["logdens"] = ld; out["z"] = zdg; out["W"] = W; out["X_p"] = xp; out["X_s"] = xs; out["X_o"] = xo; out["X_p_var"] = xpv; out["X_s_var"] = xsv; out["X_v_var"] = xvv; out["M_p"] = mp; out["M_s"] = ms;
  if (smooth) {
    std::vector<V7> Xs = Xf; std::vector<M7> Ps = Pf;
    for (int k = N - 2; k >= 0; --k) { M7 G = Pf[k] * Fk[k + 1].transpose() * pinv7(Ppr[k + 1]); Xs[k] = Xf[k] + G * (Xs[k + 1] - Xpr[k + 1]); Ps[k] = Pf[k] + G * (Ps[k + 1] - Ppr[k + 1]) * G.transpose(); }
    NumericVector sp(N), ss(N), so(N), spv(N), ssv(N), svv(N);
    for (int k = 0; k < N; ++k) { const V7& x = Xs[k]; const M7& p = Ps[k]; sp[k] = x(0) + x(5); ss[k] = x(1) + x(6); so[k] = x(2);
      spv[k] = p(0, 0) + p(5, 5) + 2.0 * p(0, 5); ssv[k] = p(1, 1) + p(6, 6) + 2.0 * p(1, 6); svv[k] = spv[k] + p(2, 2) + 2.0 * (p(0, 2) + p(2, 5)); }
    out["Xs_p"] = sp; out["Xs_s"] = ss; out["Xs_o"] = so; out["Xs_p_var"] = spv; out["Xs_s_var"] = ssv; out["Xs_v_var"] = svv;
  }
  return out;
}
