// RAE 3.0 core -- per-sample multi-component model: noise, normalisation, parameterisation,
// soft discriminators and their derivatives.  Spec steps 1-4 and 7 (per-segment part).
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#include "density.hpp"
#include "rng.hpp"
#include "types.hpp"

namespace rae {

// ---- Step 4: soft discriminators (Eq. 2) -------------------------------------------------------
// [STATED: M&M Eq. 2]  S_k(a) = 1 / (exp(-beta (a - E)) + 1)   for k != A1
inline double sigmoid(double a, double E, double beta) { return 1.0 / (std::exp(-beta * (a - E)) + 1.0); }
// [STATED: M&M Eq. 2]  S_A1(a) = Theta(a - E) * (2 / (exp(-beta (a - E)) + 1) - 1)
inline double sigmoid_a1(double a, double E, double beta) {
  return a > E ? 2.0 / (std::exp(-beta * (a - E)) + 1.0) - 1.0 : 0.0;
}
// a1_logistic: A1 as an ordinary sigmoid (preset "sequencing")  [INFERRED G-13f]
inline double score_component(int k, double a, const double* E, const double* B, bool a1_logistic = false) {
  return k == A1 && !a1_logistic ? sigmoid_a1(a, E[k], B[k]) : sigmoid(a, E[k], B[k]);
}
// ---- Step 7: derivative used for error propagation ---------------------------------------------
// [STATED: S1-C, as rendered]  dS = |beta| e^{beta(a-E)} / (e^{beta(a-E)} + 1)^2   [INFERRED G-10a, G-10d]
inline double dsigmoid(double a, double E, double beta) {
  double z = beta * (a - E);
  if (z > 700 || z < -700) return 0.0;
  double e = std::exp(z);
  return std::fabs(beta) * e / ((e + 1) * (e + 1));
}
// [STATED: S1-C]  dS_A1 = 2 Theta(a - E_A1) dS
inline double dscore_component(int k, double a, const double* E, const double* B, bool a1_logistic = false) {
  if (k == A1 && !a1_logistic) return a > E[k] ? 2.0 * dsigmoid(a, E[k], B[k]) : 0.0;
  return dsigmoid(a, E[k], B[k]);
}

// ---- G-01a: per-marker SD from the segmentation when the input has none ------------------------
// Q10 of |delta mean| * sqrt(n1 n2 / (n1 + n2)) over adjacent segment pairs, divided by c.
// Returns NaN when the sample has fewer than 10 breakpoints.  [INFERRED G-01a]
inline double sd_from_breakpoints(const Sample& s, double c) {
  std::vector<double> t;
  for (size_t i = 1; i < s.segs.size(); ++i) {
    const Segment &a = s.segs[i - 1], &b = s.segs[i];
    if (a.chrom != b.chrom || a.nmark < 1 || b.nmark < 1) continue;
    double n1 = a.nmark, n2 = b.nmark;
    t.push_back(std::fabs(a.mean - b.mean) * std::sqrt(n1 * n2 / (n1 + n2)));
  }
  if (t.size() < 10) return NaN;
  return quantile7(t, 0.10) / c;
}

inline double weighted_median_sd(const Sample& s) {
  std::vector<std::pair<double, double>> v;
  double tot = 0;
  for (const Segment& g : s.segs)
    if (!is_nan(g.sd) && g.sd > 0 && g.nmark > 1) { v.push_back({g.sd, (double)g.nmark}); tot += g.nmark; }
  if (v.empty()) return NaN;
  std::sort(v.begin(), v.end());
  double acc = 0;
  for (auto& p : v) { acc += p.second; if (acc >= tot / 2) return p.first; }
  return v.back().first;
}

// Fill in missing SDs.  `cohort_sd` is used when a sample has too few breakpoints (pass NaN on the
// first pass to learn which samples need it).
inline void resolve_sd(Sample& s, const Params& P, double cohort_sd) {
  double wm = weighted_median_sd(s);
  if (!is_nan(P.marker_sd)) { s.sd_est = P.marker_sd; }
  else if (!is_nan(wm)) { s.sd_est = wm; }
  else {
    s.flags |= F_SD_FROM_BREAKPOINTS;
    s.sd_est = sd_from_breakpoints(s, P.sd_from_breakpoints_c);
    if (is_nan(s.sd_est)) { s.sd_est = cohort_sd; s.flags |= F_SD_COHORT_MEDIAN; }
  }
}

// Local maxima over an odd window of `span` grid points: y[i] strictly above every other point in
// the window (R idiom `peaksign`, used by the original code with span 13).
inline std::vector<int> local_maxima(const std::vector<double>& y, int span) {
  std::vector<int> out;
  const int h = std::max(1, span / 2), n = (int)y.size();
  for (int i = h; i < n - h; ++i) {
    bool pk = true;
    for (int k = -h; k <= h && pk; ++k) if (k != 0 && !(y[i] > y[i + k])) pk = false;
    if (pk) out.push_back(i);
  }
  return out;
}

// [CODE D-07] Mode offset as the original code finds it: density of the marker-level segment means
// (no noise added), bandwidth DN/8, evaluated on [-DN/2, DN/2]; largest local maximum (span 11), or
// the maximum if there is none. Assumes the profile is already roughly centred.
inline double window_mode(const Sample& s, const Params& P) {
  std::vector<double> x;
  for (const Segment& g : s.segs) x.insert(x.end(), std::max(1, g.nmark), g.mean);
  const double bw = s.dn / 8.0;
  Density d = density_r(x, bw, -s.dn / 2, s.dn / 2, 512);
  std::vector<int> pk = local_maxima(d.y, 11);
  int best = (int)(std::max_element(d.y.begin(), d.y.end()) - d.y.begin());
  if (!pk.empty()) { best = pk[0]; for (int i : pk) if (d.y[i] > d.y[best]) best = i; }
  return d.x[best];
}

// ---- Preset "sequencing": detectors placed on the copy-number ladder  [INFERRED G-13] ------------
// Segment means (>= ubp_min_markers markers) are taken to linear ratio and rescaled to a marker-weighted
// mean of 1. Copy state c then sits at r(c) = r0 + c * step with step = p / (p*psi + 2(1-p)) and
// r0 = 2(1-p) / (p*psi + 2(1-p)) (purity p, mean ploidy psi).
//  - With user purity and ploidy, the ladder follows from them.
//  - Otherwise the step is the period of the segment means: for each candidate step D the comb score
//      S(D) = | sum_i w_i exp(2 pi i r_i / D) | / sum_i w_i
//    is 1 for a perfect ladder and decays with spread; harmonics D/2, D/3 ... decay faster. The largest D
//    with S >= comb_rel * max S is taken (fewest states). If that D is the top of the range (S still
//    rising: no resolvable ladder), detectors stay as RAE 2008 places them (F_STEP_FALLBACK).  [G-13b]
//    The ladder level nearest ratio 1 is the base; its copy number is fixed by taking the lowest level
//    that holds >= comb_tooth_min_frac of the genome as one copy, capped so purity <= 1.  [G-13c]
// Detectors, on the sample's centred log2 scale: A0 / D0 at round(psi) +/- 0.5 copies; A1 at the
// amplification copy number (5, or 9 above ploidy 2.7) [G-13d]; D1 at homdel_cn copies [G-13e].
// Each slope puts S = sigmoid_at_half at the neighbouring whole copy number: the paper's construction
// (S(H) = 0.001) with the ladder in place of the half-maximum point; no noise model enters.
// [INFERRED G-14a] Prior density of a mean ploidy: (1 - w) N(psi; mu1, sd) + w N(psi; mu2, sd).
inline double ploidy_prior_density(double psi, const Params& P) {
  const double k = 0.3989422804014327 / P.ploidy_prior_sd;
  auto n = [&](double mu) { const double z = (psi - mu) / P.ploidy_prior_sd; return k * std::exp(-0.5 * z * z); };
  return (1 - P.ploidy_prior_wgd) * n(P.ploidy_prior_mu1) + P.ploidy_prior_wgd * n(P.ploidy_prior_mu2);
}

inline void place_copy_step(Sample& s, const Params& P) {
  std::vector<double> r, w;
  double sw = 0, swr = 0;
  for (const Segment& g : s.segs)
    if (g.nmark >= P.ubp_min_markers) {
      double v = std::pow(2.0, g.mean);
      r.push_back(v); w.push_back((double)g.nmark); sw += g.nmark; swr += g.nmark * v;
    }
  if (r.size() < 3 || !(sw > 0)) { s.flags |= F_STEP_FALLBACK; return; }
  const double m = swr / sw;
  for (double& v : r) v /= m;
  double step, r0, purity, ploidy;
  if (!is_nan(s.in_purity) && !is_nan(s.in_ploidy) && s.in_purity > 0 && s.in_purity <= 1 && s.in_ploidy > 0) {
    purity = s.in_purity; ploidy = s.in_ploidy;
    const double den = purity * ploidy + 2 * (1 - purity);
    step = purity / den; r0 = 2 * (1 - purity) / den;
    s.step_src = 2;
  } else {
    const double TWO_PI = 6.283185307179586;
    const int nD = (int)std::floor((P.comb_step_max - P.comb_step_min) / P.comb_step_by + 0.5) + 1;
    if (nD < 2) { s.flags |= F_STEP_FALLBACK; return; }
    std::vector<double> S(nD), RE(nD), IM(nD);
    double smax = 0;
    for (int i = 0; i < nD; ++i) {
      const double D = P.comb_step_min + i * P.comb_step_by;
      double re = 0, im = 0;
      for (size_t j = 0; j < r.size(); ++j) { re += w[j] * std::cos(TWO_PI * r[j] / D); im += w[j] * std::sin(TWO_PI * r[j] / D); }
      RE[i] = re; IM[i] = im; S[i] = std::sqrt(re * re + im * im) / sw;
      if (S[i] > smax) smax = S[i];
    }
    int k = nD - 1;
    while (k > 0 && S[k] < P.comb_rel * smax) --k;
    if (k == nD - 1) { s.flags |= F_STEP_FALLBACK; return; }   // still rising at the edge
    const double D = P.comb_step_min + k * P.comb_step_by;
    const double ph = std::atan2(IM[k], RE[k]) / TWO_PI * D;   // a ladder level sits at ph (mod D)
    double t = std::fmod(ph - 1 + D / 2, D); if (t < 0) t += D;
    const double base = 1 + t - D / 2;                           // ladder level nearest ratio 1
    std::map<long, double> occ;                                  // genome fraction per level, relative to base
    for (size_t j = 0; j < r.size(); ++j) occ[std::lround((r[j] - base) / D)] += w[j] / sw;
    long jlow = 0;
    for (const auto& kv : occ) if (kv.second >= P.comb_tooth_min_frac) { jlow = kv.first; break; }
    const long kA = (long)std::floor(base / D);                  // largest admissible level (r0 >= 0)
    long kb = std::min(kA, 1 - jlow);
    if (kb < 1) kb = std::min<long>(1, kA);
    if (kb < 1) { s.flags |= F_STEP_FALLBACK; return; }
    if (P.ploidy_prior && kA > kb) {   // [INFERRED G-14a] soft prior on which level is the base: kb_hard .. kA
      const long kb_hard = kb;
      std::vector<std::pair<double, long>> cand;   // (prior density of the ploidy, base copy number)
      for (long k2 = kb_hard; k2 <= kA; ++k2) cand.push_back({ploidy_prior_density((1 - base) / D + k2, P), k2});
      std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
      kb = cand[0].second;
      s.ploidy_hard = (1 - base) / D + kb_hard;
      s.ploidy_alt = (1 - base) / D + cand[1].second;
      // The ratio is NA unless finite (runner-up underflowed to 0 or a denormal)  [G-14a]. Best density 0 too:
      // the prior gave no information, so the hard pick stands and the choice is flagged  [G-14b]
      const double q = cand[0].first / cand[1].first;
      s.ploidy_prior_ratio = std::isfinite(q) ? q : NaN;
      if (!(cand[0].first > 0) || s.ploidy_prior_ratio < P.ploidy_prior_ambig) s.flags |= F_PLOIDY_AMBIGUOUS;
    }
    step = D; r0 = base - kb * D;
    ploidy = (1 - base) / D + kb;
    purity = 2 * D / (r0 + 2 * D);
    s.comb_score = S[k];
    s.step_src = 1;
  }
  s.step = step; s.purity = purity; s.ploidy = ploidy; s.lad_m = m; s.lad_r0 = r0;
  auto x = [&](double c) { return std::log2(m * std::max(r0 + c * step, 1e-6)); };  // centred log2 of copy c
  const double L = std::log(1.0 / P.sigmoid_at_half - 1.0);
  const double kmode = std::max(1.0, std::round(ploidy));
  double amp = ploidy <= P.amp_ploidy_split ? P.amp_cn_low : P.amp_cn_high;
  if (amp < kmode + 1) amp = kmode + 1;
  s.E[A0] = x(kmode + 0.5);    s.B[A0] = L / (s.E[A0] - x(kmode));
  s.E[D0] = x(kmode - 0.5);    s.B[D0] = L / (s.E[D0] - x(kmode));
  s.E[A1] = x(amp);            s.B[A1] = L / (s.E[A1] - x(amp - 1));
  s.E[D1] = x(P.homdel_cn);    s.B[D1] = L / (s.E[D1] - x(P.homdel_cn + 0.5));
  s.a1_logistic = P.a1_logistic;
  s.flags &= ~(uint32_t)(F_A0_FALLBACK | F_D0_FALLBACK | F_D1_PCT01 | F_D1_BOUND | F_A1_BOUND);
}

// ---- Steps 1-3: normalise and parameterise one sample ------------------------------------------
// Requires s.sd_est to be set. Mutates segment means (mode-centring).
inline void parameterize_sample(Sample& s, const Params& P, uint64_t stream) {
  // "Noise added back": every probe gets a draw from Normal(segment mean, s_i), where s_i represents
  // "the intrinsic noise of the segmentation (not the probe-level noise)".  [STATED: S1-B]
  // s_i = segment SD / sqrt(markers), the standard error of the segment mean.  [INFERRED G-02g]
  // (The literal reading, s_i = SD of the probes in the segment, would give a diploid peak as wide as
  // the probe noise; Fig. S5A shows half-maxima at about +/-0.04 for probe noise of about 0.2.)
  Rng rng(P.seed, stream);  // [INFERRED G-02b]
  std::vector<double> r;
  size_t total = 0;
  for (const Segment& g : s.segs) total += std::max(1, g.nmark);
  r.reserve(total);
  for (const Segment& g : s.segs) {
    double sd = (!is_nan(g.sd) && g.sd > 0) ? g.sd : s.sd_est;
    int n = std::max(1, g.nmark);
    if (P.density_noise_sem) sd /= std::sqrt((double)n);
    for (int i = 0; i < n; ++i) r.push_back(g.mean + sd * rng.norm());
  }
  // DN  [STATED: Eq. 1] when supplied; otherwise from SD.  [INFERRED G-01b]
  if (is_nan(s.dn)) { s.dn = 0.9538726 * s.sd_est; s.flags |= F_DN_FROM_SD; }
  // Density: R density(), from -1 to 1, n = 2048, adjust = 2.  [STATED: S1-B] [INFERRED G-02a]
  double bw = !is_nan(P.density_bw) ? P.density_bw : P.density_adjust * bw_nrd0(r);
  Density d;
  int im;
  if (P.centre_window) {
    s.shift = window_mode(s, P);  // [CODE D-07]
  } else {
    // Diploid peak = global maximum; subtract its position.  [STATED: M&M] [INFERRED G-02c]
    d = density_r(r, bw, P.density_from, P.density_to, P.density_n);
    im = (int)(std::max_element(d.y.begin(), d.y.end()) - d.y.begin());
    s.shift = d.x[im];
  }
  for (Segment& g : s.segs) g.mean -= s.shift;
  for (double& v : r) v -= s.shift;
  // Density of the centred data, so that half-maxima are relative to a peak at ~0. [INFERRED G-02e]
  d = density_r(r, bw, P.density_from, P.density_to, P.density_n);
  im = (int)(std::max_element(d.y.begin(), d.y.end()) - d.y.begin());
  const std::vector<int> maxima = local_maxima(d.y, P.peak_span);  // [INFERRED G-03e] / [CODE D-05]
  if (P.diploid_nearest_zero && !maxima.empty()) {  // [CODE D-06] peak nearest the centre of the grid
    const double mid = 0.5 * (P.density_from + P.density_to);
    im = maxima[0];
    for (int i : maxima) if (std::fabs(d.x[i] - mid) < std::fabs(d.x[im] - mid)) im = i;
  }
  s.bw = bw;
  s.dens_x = d.x; s.dens_y = d.y;
  const int n = (int)d.x.size();
  const double ymax = d.y[im], xmode = d.x[im];

  // Half-maximum points: first crossing walking outward from the mode.  [STATED: S1-B] [INFERRED G-03d]
  auto halfmax = [&](int dir) {
    for (int i = im + dir; i >= 0 && i < n; i += dir) {
      if (d.y[i] < 0.5 * ymax) {
        int j = i - dir;  // last point at or above half
        double f = (d.y[j] - 0.5 * ymax) / (d.y[j] - d.y[i]);
        return d.x[j] + f * (d.x[i] - d.x[j]);
      }
    }
    s.flags |= F_NO_HALFMAX;  // [INFERRED G-03f]
    return dir > 0 ? P.density_to : P.density_from;
  };
  s.Ha = halfmax(+1);
  s.Hd = halfmax(-1);
  if (s.Ha < P.halfmax_min) s.Ha = P.halfmax_min;      // [CODE D-04]
  if (s.Hd > -P.halfmax_min) s.Hd = -P.halfmax_min;

  // Gain / loss peak: largest local maximum with 0.5 DN <= |x| <= 2 DN lying beyond the
  // half-maximum point.  [STATED: S1-B] reading [INFERRED G-03a, G-03b, G-03e]
  auto find_peak = [&](int dir, double H) {
    double best_y = -1, best_x = NaN;
    for (int i : maxima) {
      if (dir > 0 ? i <= im : i >= im) continue;
      double x = (d.x[i] - xmode) * dir;  // distance from the diploid peak on this side
      if (x < P.peak_window_lo * s.dn || x > P.peak_window_hi * s.dn) continue;
      if (x <= std::fabs(H - xmode)) continue;
      if (d.y[i] > best_y) { best_y = d.y[i]; best_x = d.x[i]; }
    }
    return best_x;
  };
  const double L = std::log(1.0 / P.sigmoid_at_half - 1.0);  // ln(999) for 1e-3
  auto set_single = [&](int k, int dir, double H, double beta_fb, uint32_t flag, double& Xout) {
    double X = find_peak(dir, H);
    Xout = X;
    if (!is_nan(X)) {
      // [STATED: S1-B] "E ... is the average of these two values" (peak and half-maximum point).
      // [CODE D-01] the code averages the loss peak with zero instead.
      s.E[k] = (k == D0 && P.d0_midpoint_zero) ? 0.5 * X : 0.5 * (X + H);
      if (P.e_inflate)  // [CODE D-08]
        for (int it = 0; it < 200 && sigmoid(H, s.E[k], dir * 100.0) > P.sigmoid_at_half; ++it) s.E[k] *= 1.1;
      s.B[k] = L / (s.E[k] - H);         // [STATED: S1-B] S(H; E, b) = 0.001
      if (P.beta_min > 0 && (std::fabs(s.B[k]) < P.beta_min || !std::isfinite(s.B[k]))) s.B[k] = dir * P.beta_min;  // [CODE D-03]
    } else {
      s.B[k] = dir * beta_fb;            // [STATED: S1-B] "fix b at a constant value of 20"
      s.E[k] = H + L / s.B[k];           // [STATED: S1-B] S(H; E, 20) = 0.001
      s.flags |= flag;
    }
  };
  set_single(A0, +1, s.Ha, P.beta_fallback_gain, F_A0_FALLBACK, s.Xa);
  set_single(D0, -1, s.Hd, P.beta_fallback_loss, F_D0_FALLBACK, s.Xd);

  // D1  [STATED: S1-B]; percentiles of mode-centred probe-level signal  [INFERRED G-04b]
  double q025, q01;
  if (!is_nan(s.q025) && !is_nan(s.q01)) { q025 = s.q025 - s.shift; q01 = s.q01 - s.shift; }
  else {  // [INFERRED G-01c] pseudo-probes stand in for probe-level data
    std::vector<double> pp;
    pp.reserve(r.size());
    for (const Segment& g : s.segs) {
      double sd = (!is_nan(g.sd) && g.sd > 0) ? g.sd : s.sd_est;
      int n = std::max(1, g.nmark);
      for (int i = 0; i < n; ++i) pp.push_back(g.mean + sd * rng.norm());  // means are already centred
    }
    std::sort(pp.begin(), pp.end());
    q025 = quantile7_sorted(pp, P.d1_pct_first);
    q01 = quantile7_sorted(pp, P.d1_pct_second);
    s.flags |= F_Q_FROM_PSEUDO;
  }
  s.q025 = q025; s.q01 = q01;  // stored on the centred scale from here on
  s.B[D1] = s.B[D0];                      // [STATED: S1-B]
  const double bound = P.d1_bound * s.E[D0];
  s.E[D1] = q025;
  if (s.E[D1] > bound) {                  // not negative enough  [INFERRED G-04a]
    s.E[D1] = q01; s.flags |= F_D1_PCT01;
    if (s.E[D1] > bound) { s.E[D1] = bound; s.flags |= F_D1_BOUND; }
  }
  // A1  [STATED: S1-B]
  s.B[A1] = P.beta_a1;
  s.E[A1] = 2 * s.E[A0] - s.Ha;           // [INFERRED G-05c] same formula after an A0 fallback
  if (s.E[A1] < P.a1_bound * s.E[A0]) { s.E[A1] = P.a1_bound * s.E[A0]; s.flags |= F_A1_BOUND; }
  if (P.copy_step) place_copy_step(s, P);   // [INFERRED G-13]
}

// Per-segment scores and squared error terms for one sample.
struct SegScores { std::vector<float> s[4], e2[4]; };  // e2 = (dS * sd / sqrt(n))^2
inline SegScores score_segments(const Sample& smp) {
  SegScores out;
  size_t n = smp.segs.size();
  for (int k = 0; k < 4; ++k) { out.s[k].resize(n); out.e2[k].resize(n); }
  for (size_t i = 0; i < n; ++i) {
    const Segment& g = smp.segs[i];
    double sd = (!is_nan(g.sd) && g.sd > 0) ? g.sd : smp.sd_est;
    double se = sd / std::sqrt((double)std::max(1, g.nmark));
    for (int k = 0; k < 4; ++k) {
      out.s[k][i] = (float)score_component(k, g.mean, smp.E, smp.B, smp.a1_logistic);
      double t = dscore_component(k, g.mean, smp.E, smp.B, smp.a1_logistic) * se;
      out.e2[k][i] = (float)(t * t);
    }
  }
  return out;
}

}  // namespace rae
