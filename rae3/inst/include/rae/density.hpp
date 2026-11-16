// RAE 3.0 core -- kernel density estimate reproducing R's density() and small statistics helpers.
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "types.hpp"

namespace rae {

// Quantile, R type 7 (the R default).  [INFERRED G-04c]  `v` is sorted in place.
inline double quantile7(std::vector<double>& v, double p) {
  if (v.empty()) return NaN;
  std::sort(v.begin(), v.end());
  double h = (v.size() - 1) * p;
  size_t lo = (size_t)std::floor(h);
  size_t hi = std::min(lo + 1, v.size() - 1);
  return v[lo] + (h - lo) * (v[hi] - v[lo]);
}
inline double quantile7_sorted(const std::vector<double>& v, double p) {
  if (v.empty()) return NaN;
  double h = (v.size() - 1) * p;
  size_t lo = (size_t)std::floor(h);
  size_t hi = std::min(lo + 1, v.size() - 1);
  return v[lo] + (h - lo) * (v[hi] - v[lo]);
}
inline double median(std::vector<double> v) { return quantile7(v, 0.5); }

// R's bw.nrd0: 0.9 * min(sd, IQR/1.34) * n^(-1/5)
inline double bw_nrd0(const std::vector<double>& x) {
  size_t n = x.size();
  if (n < 2) return 1.0;
  double m = 0;
  for (double v : x) m += v;
  m /= n;
  double ss = 0;
  for (double v : x) ss += (v - m) * (v - m);
  double hi = std::sqrt(ss / (n - 1));
  std::vector<double> s(x);
  std::sort(s.begin(), s.end());
  double iqr = quantile7_sorted(s, 0.75) - quantile7_sorted(s, 0.25);
  double lo = std::min(hi, iqr / 1.34);
  if (lo == 0) { lo = hi; if (lo == 0) lo = std::fabs(x[0]); if (lo == 0) lo = 1.0; }
  return 0.9 * lo * std::pow((double)n, -0.2);
}

struct Density { std::vector<double> x, y; double bw = NaN; };

// Gaussian KDE on n equally spaced points over [from, to], following stats::density.default:
// linear binning of the data on [from - 4bw, to + 4bw], convolution with the kernel, then linear
// interpolation onto the output grid. R does the convolution by FFT; here it is done directly with
// the kernel truncated at 6 bw, which agrees with R to well below plotting precision. [INFERRED G-02d]
// Points outside the binning range still count in the normalising total, as in R. [INFERRED G-02f]
inline Density density_r(const std::vector<double>& data, double bw, double from, double to, int n) {
  Density d;
  d.bw = bw;
  const double lo = from - 4 * bw, up = to + 4 * bw;
  const double xdelta = (up - lo) / (n - 1);
  std::vector<double> y(2 * n, 0.0);
  const double w = 1.0 / (double)data.size();
  const int ixmax = n - 2;
  for (double v : data) {
    double xpos = (v - lo) / xdelta;
    if (!(xpos > -2 && xpos < n + 2)) continue;
    int ix = (int)std::floor(xpos);
    double fx = xpos - ix;
    if (ix >= 0 && ix <= ixmax) { y[ix] += w * (1 - fx); y[ix + 1] += w * fx; }
    else if (ix == -1) y[0] += w * fx;
    else if (ix == ixmax + 1) y[ix] += w * (1 - fx);
  }
  int half = (int)std::ceil(6 * bw / xdelta);
  std::vector<double> kern(half + 1);
  const double c = 1.0 / (bw * std::sqrt(2 * M_PI));
  for (int k = 0; k <= half; ++k) { double z = k * xdelta / bw; kern[k] = c * std::exp(-0.5 * z * z); }
  std::vector<double> dens(n, 0.0);
  for (int h = 0; h < n; ++h) {
    if (y[h] == 0) continue;
    int g0 = std::max(0, h - half), g1 = std::min(n - 1, h + half);
    for (int g = g0; g <= g1; ++g) dens[g] += y[h] * kern[std::abs(g - h)];
  }
  d.x.resize(n); d.y.resize(n);
  for (int i = 0; i < n; ++i) {
    double xo = from + (to - from) * i / (n - 1);
    double pos = (xo - lo) / xdelta;
    int ix = std::min(std::max((int)std::floor(pos), 0), n - 2);
    double f = pos - ix;
    d.x[i] = xo;
    d.y[i] = dens[ix] * (1 - f) + dens[ix + 1] * f;
  }
  return d;
}

// Benjamini-Hochberg q-values (identical to R p.adjust(method = "BH")).  [STATED: Main]
inline std::vector<double> bh_adjust(const std::vector<double>& p) {
  size_t n = p.size();
  std::vector<size_t> o(n);
  for (size_t i = 0; i < n; ++i) o[i] = i;
  std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return p[a] > p[b]; });
  std::vector<double> q(n);
  double running = 1.0;
  for (size_t r = 0; r < n; ++r) {  // r = 0 is the largest p (rank n)
    double v = p[o[r]] * (double)n / (double)(n - r);
    running = std::min(running, v);
    q[o[r]] = running;
  }
  return q;
}

}  // namespace rae
