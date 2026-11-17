// experiments/features/features.cpp -- sample entropy and permutation entropy (EXPERIMENTAL; compiled by lib.R with Rcpp::sourceCpp)
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
using namespace Rcpp;

// Sample entropy SampEn(m, r) = -ln(A / B): B = pairs of templates of length m within r (Chebyshev), A = the same for length m + 1; self-matches excluded,
// the same N - m templates for both lengths (Richman & Moorman 2000). NA when A or B is 0.
// [[Rcpp::export]]
double sampen_cpp(NumericVector x, int m, double r) {
  int n = x.size(); if (n <= m + 1) return NA_REAL;
  double A = 0, B = 0; int nt = n - m;
  for (int i = 0; i < nt - 1; i++) {
    for (int j = i + 1; j < nt; j++) {
      bool ok = true;
      for (int k = 0; k < m; k++) if (std::fabs(x[i + k] - x[j + k]) > r) { ok = false; break; }
      if (!ok) continue;
      B += 1;
      if (std::fabs(x[i + m] - x[j + m]) <= r) A += 1;
    }
  }
  if (A == 0 || B == 0) return NA_REAL;
  return -std::log(A / B);
}

// Normalized permutation entropy (Bandt & Pompe 2002) of order `ord` and delay `tau`: Shannon entropy of the ordinal patterns / ln(ord!).
// [[Rcpp::export]]
double permen_cpp(NumericVector x, int ord, int tau) {
  int n = x.size(); int nv = n - (ord - 1) * tau; if (nv < 2) return NA_REAL;
  std::map<long, double> counts; std::vector<int> idx(ord); std::vector<double> v(ord);
  for (int i = 0; i < nv; i++) {
    for (int k = 0; k < ord; k++) { v[k] = x[i + k * tau]; idx[k] = k; }
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return v[a] < v[b]; });
    long code = 0; for (int k = 0; k < ord; k++) code = code * ord + idx[k];
    counts[code] += 1;
  }
  double h = 0; for (auto &c : counts) { double p = c.second / nv; h -= p * std::log(p); }
  double f = 1; for (int k = 2; k <= ord; k++) f *= k;
  return h / std::log(f);
}
