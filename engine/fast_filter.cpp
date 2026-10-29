// Fast port of run_gls_filter() from model_engine.R (same recursion, same
// guards). Used only to speed up multi-fit analyses; it is verified against the
// R engine in tests/test_pipeline.R and never replaces it for reported fits.
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

struct OUBlock { double phi_int, q11, q13, q33; };

// Exact OU integrals over dt (mirrors .compute_ou_block, incl. Taylor branch).
static inline OUBlock ou_block(double k, double dt, double V) {
  double x = k * dt, s2 = 2.0 * k * V;
  OUBlock b;
  if (x < 1e-3) {
    b.phi_int = dt - (k * dt * dt) / 2.0 + (k * k * dt * dt * dt) / 6.0;
    b.q11 = s2 * (dt - k * dt * dt + (2.0 * k * k * std::pow(dt, 3)) / 3.0);
    b.q13 = s2 * ((dt * dt) / 2.0 - (k * std::pow(dt, 3)) / 2.0 + (7.0 * k * k * std::pow(dt, 4)) / 24.0);
    b.q33 = s2 * (std::pow(dt, 3) / 3.0 - (k * std::pow(dt, 4)) / 4.0 + (7.0 * k * k * std::pow(dt, 5)) / 60.0);
  } else {
    double e1 = std::exp(-x), e2 = std::exp(-2.0 * x);
    b.phi_int = (1.0 - e1) / k;
    b.q11 = s2 * (1.0 - e2) / (2.0 * k);
    b.q13 = s2 * (1.0 - 2.0 * e1 + e2) / (2.0 * k * k);
    b.q33 = s2 * (2.0 * x - 3.0 + 4.0 * e1 - e2) / (2.0 * k * k * k);
  }
  return b;
}

struct PassOut {
  std::vector<double> S, Sung, v0, vv, x0, xv, pvar;  // x0/xv: N x 4 (row-major); pvar: N x 4 diag + P34
  double nu0;
};

static PassOut gls_pass(const std::vector<double>& dy, const std::vector<double>& gm,
                        double kp, double ks, double lp, double lR, bool keep) {
  const int N = dy.size();
  PassOut o;
  o.S.assign(N, 0); o.Sung.assign(N, 0); o.v0.assign(N, 0); o.vv.assign(N, 0);
  if (keep) { o.x0.assign(4 * N, 0); o.xv.assign(4 * N, 0); o.pvar.assign(5 * N, 0); }

  double X0[4] = {0, 0, 0, 0}, Xv[4] = {0, 0, 0, 0};
  double P[4][4] = {{lp, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

  for (int k = 0; k < N; ++k) {
    double dt = dy[k];
    double dyn_lR = std::min(lR * gm[k], 1e12);
    OUBlock bs = ou_block(ks, dt, 1.0);
    double phi_p = ou_block(kp, dt, 1.0).phi_int;

    double Phi[4][4] = {{0}};
    Phi[0][0] = std::exp(-kp * dt); Phi[2][0] = phi_p;
    Phi[1][1] = std::exp(-ks * dt); Phi[3][1] = bs.phi_int;

    double X0p[4], Xvp[4];
    for (int i = 0; i < 4; ++i) {
      X0p[i] = 0; Xvp[i] = 0;
      for (int j = 0; j < 4; ++j) { X0p[i] += Phi[i][j] * X0[j]; Xvp[i] += Phi[i][j] * Xv[j]; }
    }

    OUBlock bp = ou_block(kp, dt, lp);

    double Q[4][4] = {{0}};
    Q[0][0] = bp.q11; Q[0][2] = bp.q13; Q[2][0] = bp.q13; Q[2][2] = bp.q33;
    Q[1][1] = bs.q11; Q[1][3] = bs.q13; Q[3][1] = bs.q13; Q[3][3] = bs.q33;

    double T[4][4], Pp[4][4];
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
      double s = 0; for (int l = 0; l < 4; ++l) s += Phi[i][l] * P[l][j]; T[i][j] = s;
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
      double s = 0; for (int l = 0; l < 4; ++l) s += T[i][l] * Phi[j][l]; Pp[i][j] = s + Q[i][j];
    }

    double HPH = Pp[2][2] + Pp[3][3] - 2.0 * Pp[2][3];
    double S = std::max(HPH + dyn_lR * dt, 1e-12);
    double Kg[4];
    for (int i = 0; i < 4; ++i) Kg[i] = (-Pp[i][2] + Pp[i][3]) / S;
    o.Sung[k] = std::max(HPH + lR * dt, 1e-12);

    double v0 = 1.0 - (-X0p[2] + X0p[3]);
    double vv = dt - (-Xvp[2] + Xvp[3]);
    for (int i = 0; i < 4; ++i) { X0[i] = X0p[i] + Kg[i] * v0; Xv[i] = Xvp[i] + Kg[i] * vv; }

    // Joseph update: (I - K H) Pp (I - K H)' + K (dyn_lR dt) K'
    double IKH[4][4];
    const double Hh[4] = {0, 0, -1, 1};
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) IKH[i][j] = (i == j ? 1.0 : 0.0) - Kg[i] * Hh[j];
    double A1[4][4], Pn[4][4];
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
      double s = 0; for (int l = 0; l < 4; ++l) s += IKH[i][l] * Pp[l][j]; A1[i][j] = s;
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
      double s = 0; for (int l = 0; l < 4; ++l) s += A1[i][l] * IKH[j][l];
      Pn[i][j] = s + Kg[i] * (dyn_lR * dt) * Kg[j];
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) P[i][j] = 0.5 * (Pn[i][j] + Pn[j][i]);

    if (keep) {
      for (int i = 0; i < 4; ++i) { o.x0[4 * k + i] = X0[i]; o.xv[4 * k + i] = Xv[i]; o.pvar[5 * k + i] = P[i][i]; }
      o.pvar[5 * k + 4] = P[2][3];
    }
    // IPFM reset of the phase states at the beat boundary
    X0[2] = X0[3] = 0; Xv[2] = Xv[3] = 0;
    for (int j = 0; j < 4; ++j) { P[2][j] = 0; P[3][j] = 0; P[j][2] = 0; P[j][3] = 0; }

    o.S[k] = S; o.v0[k] = v0; o.vv[k] = vv;
  }
  double den = 0, num = 0;
  for (int k = 0; k < N; ++k) { den += o.vv[k] * o.vv[k] / o.S[k]; num += o.v0[k] * o.vv[k] / o.S[k]; }
  o.nu0 = (std::isfinite(den) && den > 1e-12) ? num / den : 0.0;
  return o;
}

// [[Rcpp::export]]
List gls_filter_cpp(NumericVector dy_, double kp, double ks, double lp, double lR,
                    double jump_threshold, double jump_power, double lik_gate_max, bool full) {
  const int N = dy_.size();
  std::vector<double> dy(dy_.begin(), dy_.end()), gm(N, 1.0);

  // Artifact multiplier (mirrors the R engine, incl. one-sided ends)
  for (int k = 0; k < N; ++k) {
    double jb = 0.0, jf = 0.0;
    if (k > 0) jb = std::fabs(dy[k] - dy[k - 1]) / std::max(dy[k - 1], 1e-6);
    if (k < N - 1) jf = std::fabs(dy[k + 1] - dy[k]) / std::max(dy[k], 1e-6);
    double J = std::max(jb, jf);
    double m = std::min(1.0 + std::pow(J / jump_threshold, jump_power), 1e10);
    if (!std::isfinite(m)) m = 1e10;
    gm[k] = m;
  }

  PassOut res = gls_pass(dy, gm, kp, ks, lp, lR, full);

  if (!std::isfinite(res.nu0) || res.nu0 <= 0)
    return List::create(_["LL_star"] = R_NegInf, _["nu0"] = NA_REAL, _["sig2"] = NA_REAL);

  std::vector<double> v(N);
  int n_eff = 0; for (int k = 0; k < N; ++k) if (gm[k] < lik_gate_max) ++n_eff;
  bool all = n_eff < 2; if (all) n_eff = N;
  double s2 = 0, sumlogS = 0;
  for (int k = 0; k < N; ++k) {
    v[k] = res.v0[k] - res.nu0 * res.vv[k];
    if (all || gm[k] < lik_gate_max) { s2 += v[k] * v[k] / res.S[k]; sumlogS += std::log(res.S[k]); }
  }
  double sig2 = std::max(s2 / n_eff, 1e-12);
  double LL = -0.5 * sumlogS - (n_eff / 2.0) * std::log(sig2);

  List out = List::create(_["LL_star"] = LL, _["nu0"] = res.nu0, _["sig2"] = sig2);
  if (full) {
    NumericVector xp(N), xs(N), rr(N), zin(N), zdg(N), xpv(N), xsv(N);
    double den = 0; for (int k = 0; k < N; ++k) den += res.vv[k] * res.vv[k] / res.S[k];
    double inv_den = (std::isfinite(den) && den > 1e-12) ? 1.0 / den : NA_REAL;
    for (int k = 0; k < N; ++k) {
      xp[k] = res.x0[4 * k] - res.nu0 * res.xv[4 * k];
      xs[k] = res.x0[4 * k + 1] - res.nu0 * res.xv[4 * k + 1];
      double ph_p = res.x0[4 * k + 2] - res.nu0 * res.xv[4 * k + 2];
      double ph_s = res.x0[4 * k + 3] - res.nu0 * res.xv[4 * k + 3];
      rr[k] = (1.0 + ph_p - ph_s) / res.nu0;
      zin[k] = v[k] / std::sqrt(res.S[k] * sig2);
      zdg[k] = v[k] / std::sqrt(res.Sung[k] * sig2);
      xpv[k] = (res.pvar[5 * k] + res.xv[4 * k] * res.xv[4 * k] * inv_den) * sig2;
      xsv[k] = (res.pvar[5 * k + 1] + res.xv[4 * k + 1] * res.xv[4 * k + 1] * inv_den) * sig2;
    }
    out["X_p"] = xp; out["X_s"] = xs; out["X_p_var"] = xpv; out["X_s_var"] = xsv;
    out["implied_rr"] = rr; out["std_innov"] = zin; out["innovations_diagnostic"] = zdg;
    out["gate_multiplier"] = NumericVector(gm.begin(), gm.end());
  }
  return out;
}
