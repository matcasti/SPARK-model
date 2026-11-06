// EXPERIMENTAL. Regime-switching targets ("latent Markov states") for SPARK. The two OU drives relax (rates kappa_p, kappa_s) towards a target (mu_p, mu_s)
// that depends on a hidden regime s_k in {0..K-1} which is constant within each beat interval and follows a Markov chain over beats:
//   x_b(t+h) = e_b x_b(t) + (1 - e_b) mu_b(s) + OU noise,   I_net = -I_p + I_s,   y_k = 1 - nu0 dt_k = I_net + eta_k.
// Covariances do not depend on the regime, so the regimes only shift the predicted mean: the filter is a Gaussian mixture filter with collapse after every beat (GPB1),
// whose log-likelihood increment is the log of the mixture predictive density. nu0 and sigma^2 are explicit (not profiled): the mean depends on the regime path.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
typedef Eigen::Matrix<double, 3, 3> M3; typedef Eigen::Matrix<double, 3, 1> V3;
struct OU { double phi, q11, q13, q33; };
static inline OU ou(double k, double dt, double V) {
  double x = k * dt, s2 = 2.0 * k * V; OU b;
  if (x < 1e-3) { b.phi = dt - (k * dt * dt) / 2.0 + (k * k * dt * dt * dt) / 6.0; b.q11 = s2 * (dt - k * dt * dt + (2.0 * k * k * std::pow(dt, 3)) / 3.0);
    b.q13 = s2 * ((dt * dt) / 2.0 - (k * std::pow(dt, 3)) / 2.0 + (7.0 * k * k * std::pow(dt, 4)) / 24.0); b.q33 = s2 * (std::pow(dt, 3) / 3.0 - (k * std::pow(dt, 4)) / 4.0 + (7.0 * k * k * std::pow(dt, 5)) / 60.0); }
  else { double e1 = std::exp(-x), e2 = std::exp(-2.0 * x); b.phi = (1.0 - e1) / k; b.q11 = s2 * (1.0 - e2) / (2.0 * k); b.q13 = s2 * (1.0 - 2.0 * e1 + e2) / (2.0 * k * k); b.q33 = s2 * (2.0 * x - 3.0 + 4.0 * e1 - e2) / (2.0 * k * k * k); }
  return b;
}
// [[Rcpp::export]]
List reg_filter(NumericVector dy_, NumericVector gm_, double kp, double ks, double lp, double lR, double sig2, double nu0,
                NumericMatrix mu, NumericMatrix Pi, NumericVector w0, double lik_gate_max, bool full) {
  const int N = dy_.size(), K = mu.nrow();
  V3 X = V3::Zero(); M3 P = M3::Zero(); P(0, 0) = lp * sig2; P(1, 1) = sig2;
  std::vector<double> w(K); for (int j = 0; j < K; ++j) w[j] = w0[j];
  double LL = 0; NumericVector ld(N), zdg(N), xp(N), xs(N), xpv(N), xsv(N), reg_mean(N); NumericMatrix W(N, K);
  M3 I3 = M3::Identity(); std::vector<double> a(K), lj(K), wn(K); std::vector<V3> mj(K);
  for (int k = 0; k < N; ++k) {
    double dt = dy_[k], R = std::min(lR * gm_[k], 1e12) * dt * sig2, Rung = lR * dt * sig2;
    OU bp = ou(kp, dt, lp * sig2), bs = ou(ks, dt, sig2); double ep = std::exp(-kp * dt), es = std::exp(-ks * dt);
    M3 F = M3::Zero(), Q = M3::Zero(); F(0, 0) = ep; F(1, 1) = es; F(2, 0) = -bp.phi; F(2, 1) = bs.phi;
    Q(0, 0) = bp.q11; Q(0, 2) = Q(2, 0) = -bp.q13; Q(1, 1) = bs.q11; Q(1, 2) = Q(2, 1) = bs.q13; Q(2, 2) = bp.q33 + bs.q33;
    M3 Pp = F * P * F.transpose() + Q; double S = std::max(Pp(2, 2) + R, 1e-12), Sung = std::max(Pp(2, 2) + Rung, 1e-12); V3 Kg = Pp.col(2) / S; double y = 1.0 - nu0 * dt;
    V3 Xp = F * X; double ymean = 0, mx = -1e300;
    for (int j = 0; j < K; ++j) { double s = 0; for (int i = 0; i < K; ++i) s += w[i] * Pi(i, j); a[j] = s;
      V3 b; b << (1.0 - ep) * mu(j, 0), (1.0 - es) * mu(j, 1), -(dt - bp.phi) * mu(j, 0) + (dt - bs.phi) * mu(j, 1); mj[j] = Xp + b; ymean += a[j] * mj[j](2);
      double v = y - mj[j](2); lj[j] = std::log(std::max(a[j], 1e-300)) - 0.5 * v * v / S - 0.5 * std::log(2.0 * M_PI * S); if (lj[j] > mx) mx = lj[j]; }
    double sum = 0; for (int j = 0; j < K; ++j) sum += std::exp(lj[j] - mx); double logp = mx + std::log(sum); ld[k] = logp; if (gm_[k] < lik_gate_max) LL += logp;
    for (int j = 0; j < K; ++j) wn[j] = std::exp(lj[j] - logp);
    double spread_pred = 0; for (int j = 0; j < K; ++j) spread_pred += a[j] * (mj[j](2) - ymean) * (mj[j](2) - ymean); zdg[k] = (y - ymean) / std::sqrt(Sung + spread_pred);
    V3 Xn = V3::Zero(); std::vector<V3> mm(K); for (int j = 0; j < K; ++j) { mm[j] = mj[j] + Kg * (y - mj[j](2)); Xn += wn[j] * mm[j]; }
    M3 IKH = I3; IKH.col(2) -= Kg; M3 Pn = IKH * Pp * IKH.transpose() + Kg * Kg.transpose() * R; for (int j = 0; j < K; ++j) { V3 d = mm[j] - Xn; Pn += wn[j] * d * d.transpose(); }
    X = Xn; P = 0.5 * (Pn + Pn.transpose()); for (int j = 0; j < K; ++j) { w[j] = wn[j]; W(k, j) = wn[j]; }
    xp[k] = X(0); xs[k] = X(1); xpv[k] = P(0, 0); xsv[k] = P(1, 1);
  }
  List out = List::create(_["LL"] = LL); if (!full) return out;
  out["logdens"] = ld; out["z"] = zdg; out["X_p"] = xp; out["X_s"] = xs; out["X_p_var"] = xpv; out["X_s_var"] = xsv; out["W"] = W; return out;
}
