// RAE 3.0 core -- chromosome-arm geometry, the split of segments at the centromere, and per-sample arm
// baselines (the copy level an arm sits at when nothing is focally altered). Shared by the opt-in
// focal-only analysis and the arm-level test (arm_test, below). Nothing here is called unless one of those is switched on.
//
// The two-level fit follows the description of the broad-level ("ziggurat") deconvolution in Mermel CH et al.
// (2011) Genome Biol 12:R41, Supplementary Methods: each chromosome is described by one level or two, chosen
// by a BIC with penalty (2k - 1) ln n. Reimplemented from that description alone.  [INFERRED G-14c]
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include "density.hpp"
#include "model.hpp"
#include "rng.hpp"
#include "types.hpp"

namespace rae {

struct Arm {
  int chrom;                   // 0-based chromosome index
  char arm;                    // 'p', 'q', or 'w' (whole chromosome: centromere unknown)
  int32_t m_first, m_last;     // global marker indices, inclusive
  int64_t bp_start, bp_end;    // from the first and last marker
  bool eligible;               // at least arm_min_markers markers
};

// Arms of every chromosome that has markers, ordered by marker index. A marker belongs to p when its midpoint
// is left of the centromere start, to q when right of its end; inside the centromere it goes to the nearer
// side (a tie goes to q). An arm with no markers is not listed. Unknown centromere: one arm 'w'.
// [INFERRED G-14c]
inline std::vector<Arm> arm_layout(const Genome& G, const Params& P) {
  std::vector<Arm> out;
  for (int c = 0; c < G.n_chrom; ++c) {
    const int lo = G.chrom_first[c], hi = G.chrom_first[c + 1];   // [lo, hi)
    if (lo >= hi) continue;
    auto add = [&](char a, int first, int last) {
      if (first > last) return;
      out.push_back({c, a, first, last, G.m_start[first], G.m_end[last], last - first + 1 >= P.arm_min_markers});
    };
    if (G.cen_start[c] <= 0 || G.cen_end[c] <= 0) { add('w', lo, hi - 1); continue; }
    const double cs = (double)G.cen_start[c], ce = (double)G.cen_end[c];
    int k = lo;   // first q marker
    while (k < hi) {
      const double mid = 0.5 * ((double)G.m_start[k] + (double)G.m_end[k]);
      const bool q = mid > ce || (mid >= cs && mid - cs >= ce - mid);
      if (q) break;
      ++k;
    }
    add('p', lo, k - 1);
    add('q', k, hi - 1);
  }
  return out;
}

// Cut every segment that crosses the p/q boundary into two, at the boundary. nmark is divided by each side's
// share of the grid markers (rounded, at least 1 each, the total kept). The segment's outer bp coordinates are
// kept as given; the cut is at one bp position, cut_bp = (start of the first q marker) - 1: the p piece ends at
// cut_bp and the q piece starts at cut_bp + 1, so the two lengths sum to the segment's. A segment whose cut_bp is
// not strictly inside its own coordinates is left whole. Every segment's `orig` is set to its input index, so a
// cut pair can be counted once (n_altered). mean and sd are copied. A segment with fewer than 2 markers cannot become two real pieces and is
// left whole; it is then judged on the arm that holds its start marker.  [INFERRED G-14c]
inline void split_segments_at_centromeres(std::vector<Sample>& S, const Genome& G, const std::vector<Arm>& arms) {
  std::vector<int32_t> cut(G.n_chrom, -1);   // last p marker of a chromosome that has both arms
  for (size_t i = 0; i + 1 < arms.size(); ++i)
    if (arms[i].chrom == arms[i + 1].chrom && arms[i].arm == 'p' && arms[i + 1].arm == 'q') cut[arms[i].chrom] = arms[i].m_last;
  for (Sample& s : S) {
    std::vector<Segment> out;
    out.reserve(s.segs.size() + 4);
    for (size_t i = 0; i < s.segs.size(); ++i) {
      Segment g = s.segs[i];
      g.orig = (int32_t)i;
      const int32_t b = g.chrom >= 0 && g.chrom < G.n_chrom ? cut[g.chrom] : -1;
      if (b < 0 || !(g.start <= b && b < g.end) || g.nmark < 2) { out.push_back(g); continue; }
      const int64_t cut_bp = G.m_start[b + 1] - 1;
      if (!(g.bp_start <= cut_bp && cut_bp < g.bp_end)) { out.push_back(g); continue; }
      Segment l = g, r = g;
      l.end = b; r.start = b + 1;
      const double nl = (double)(b - g.start + 1), tot = (double)(g.end - g.start + 1);
      l.nmark = (int32_t)std::min<long long>(std::max<long long>(std::llround(g.nmark * nl / tot), 1), g.nmark - 1);
      r.nmark = g.nmark - l.nmark;
      l.bp_end = cut_bp;
      r.bp_start = cut_bp + 1;
      out.push_back(l); out.push_back(r);
    }
    s.segs.swap(out);
  }
}

// Index of the arm holding the segment's start marker, or -1.
inline int arm_index(const std::vector<Arm>& arms, const Segment& g) {
  int a = (int)(std::upper_bound(arms.begin(), arms.end(), g.start,
                                 [](int32_t m, const Arm& x) { return m < x.m_first; }) - arms.begin()) - 1;
  if (a < 0 || arms[a].chrom != g.chrom || g.start > arms[a].m_last) a = -1;
  return a;
}

// Copy units apply when the sample has a usable ladder (G-13); otherwise the centred log2 value is used.
inline bool uses_copy_units(const Sample& s) {
  return s.step_src > 0 && !(s.flags & F_STEP_FALLBACK) && std::isfinite(s.step) && s.step > 0 &&
         std::isfinite(s.lad_m) && s.lad_m > 0 && std::isfinite(s.lad_r0);
}
// The segment's value: copy number on the ladder, c = (2^mean / lad_m - lad_r0) / step, or its log2 mean.
inline double sample_value(const Sample& s, const Segment& g, bool copy_units) {
  return copy_units ? (std::exp2(g.mean) / s.lad_m - s.lad_r0) / s.step : g.mean;
}

// Lower weighted median: the smallest value at which the cumulative weight reaches half the total. Used by both
// estimators, so a tie between two values resolves to the smaller one.  `x` is (value, weight), any order.
inline double lower_weighted_median(std::vector<std::pair<double, double>> x) {
  if (x.empty()) return NaN;
  std::sort(x.begin(), x.end());
  double tot = 0;
  for (const auto& p : x) tot += p.second;
  const double half = 0.5 * tot * (1 - 1e-12);   // guards the comparison against rounding in the cumulative sum
  double cum = 0;
  for (const auto& p : x) { cum += p.second; if (cum >= half) return p.first; }
  return x.back().first;
}

// One or two levels for a chromosome's segments (in positional order), by L1 change-point fit.
// Weights are the markers, rescaled to sum to n (the number of segments). Model 1: one level, the weighted
// median; model 2: a level either side of one breakpoint t between segments. BIC_k = 2n ln(SAD_k / n + 1e-12) +
// pen_k ln n with pen_1 = 1, pen_2 = 3, i.e. (2k - 1) ln n; the smaller wins and a tie goes to model 1.
// Mermel 2011 Supplementary Methods.  [INFERRED G-14c]
struct TwoLevel { int t; double mu[2]; };   // segments [0, t) at mu[0], [t, n) at mu[1]; t == n: one level
inline TwoLevel two_level_fit(const std::vector<double>& v, const std::vector<double>& nmark) {
  const int n = (int)v.size();
  double sw = 0;
  for (double x : nmark) sw += x;
  std::vector<double> w(n);
  for (int i = 0; i < n; ++i) w[i] = nmark[i] * n / sw;
  std::vector<int> ord(n);
  for (int i = 0; i < n; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](int a, int b) { return v[a] < v[b]; });
  // weighted median and SAD of the segments with index in [a, b), scanning in value order
  auto level = [&](int a, int b, double& mu, double& sad) {
    double tot = 0;
    for (int i = a; i < b; ++i) tot += w[i];
    const double half = 0.5 * tot * (1 - 1e-12);
    double cum = 0;
    mu = NaN;
    for (int i : ord) {
      if (i < a || i >= b) continue;
      cum += w[i];
      if (cum >= half) { mu = v[i]; break; }
    }
    sad = 0;
    for (int i = a; i < b; ++i) sad += w[i] * std::fabs(v[i] - mu);
  };
  TwoLevel best{n, {NaN, NaN}};
  double sad1;
  level(0, n, best.mu[0], sad1);
  best.mu[1] = best.mu[0];
  const double ln_n = std::log((double)n);
  double bic = 2.0 * n * std::log(sad1 / n + 1e-12) + 1.0 * ln_n;
  for (int t = 1; t < n; ++t) {
    double ml, mr, sl, sr;
    level(0, t, ml, sl);
    level(t, n, mr, sr);
    const double b2 = 2.0 * n * std::log((sl + sr) / n + 1e-12) + 3.0 * ln_n;
    if (b2 < bic) { bic = b2; best = {t, {ml, mr}}; }
  }
  return best;
}

// Per-arm levels from per-segment data. `arm` is the arm index of each segment (-1: ignore), `chrom` its
// chromosome, in positional order within a chromosome. method 0: the weighted median on the arm; 1: the two-level
// fit per chromosome, each arm taking the level that covers most of its weight. NaN where an arm has no data.
inline std::vector<double> arm_levels(const std::vector<int>& chrom, const std::vector<int>& arm, const std::vector<double>& v,
                                      const std::vector<double>& w, int n_arms, int method) {
  std::vector<double> lev(n_arms, NaN);
  const size_t n = v.size();
  auto ok = [&](size_t i) { return arm[i] >= 0 && std::isfinite(v[i]) && std::isfinite(w[i]) && w[i] > 0; };
  if (method == 0) {
    std::vector<std::vector<std::pair<double, double>>> g(n_arms);
    for (size_t i = 0; i < n; ++i) if (ok(i)) g[arm[i]].push_back({v[i], w[i]});
    for (int a = 0; a < n_arms; ++a) lev[a] = lower_weighted_median(g[a]);
    return lev;
  }
  std::vector<int> cs;   // distinct chromosomes
  for (size_t i = 0; i < n; ++i) if (ok(i)) cs.push_back(chrom[i]);
  std::sort(cs.begin(), cs.end());
  cs.erase(std::unique(cs.begin(), cs.end()), cs.end());
  std::vector<std::vector<double>> wt(n_arms, std::vector<double>(2, 0.0));
  std::vector<std::vector<double>> mu(n_arms, std::vector<double>(2, NaN));
  std::vector<char> seen(n_arms, 0);
  for (int c : cs) {
    std::vector<size_t> idx;
    std::vector<double> vv, ww;
    for (size_t i = 0; i < n; ++i) if (ok(i) && chrom[i] == c) { idx.push_back(i); vv.push_back(v[i]); ww.push_back(w[i]); }
    const TwoLevel f = two_level_fit(vv, ww);
    for (size_t j = 0; j < idx.size(); ++j) {
      const int a = arm[idx[j]], k = (int)j < f.t ? 0 : 1;
      wt[a][k] += w[idx[j]];
      mu[a][k] = f.mu[k]; seen[a] = 1;
    }
  }
  for (int a = 0; a < n_arms; ++a)
    if (seen[a]) lev[a] = wt[a][1] > wt[a][0] ? mu[a][1] : mu[a][0];   // a tie keeps the left level
  return lev;
}

struct ArmLevel { double level; double coverage; };

// Per arm: the baseline level and the share of the arm's markers whose segment sits at it. Copy units: within
// 0.5 copies of the level. Log2 scale: the residual (segment mean minus level, on the neutral-0 scale) is below 0.5
// for both the gain (A0) and the loss (D0) detector. Segments should already be split at the centromere.
inline std::vector<ArmLevel> arm_baselines(const Sample& s, const std::vector<Arm>& arms, const Params& P) {
  const int na = (int)arms.size();
  const bool cu = uses_copy_units(s);
  std::vector<int> chrom, arm;
  std::vector<double> v, w;
  for (const Segment& g : s.segs) {
    const int a = arm_index(arms, g);
    chrom.push_back(g.chrom); arm.push_back(a);
    v.push_back(sample_value(s, g, cu)); w.push_back((double)g.nmark);
  }
  const std::vector<double> lev = arm_levels(chrom, arm, v, w, na, P.arm_baseline);
  const bool have_det = cu || (std::isfinite(s.E[A0]) && std::isfinite(s.E[D0]) && std::isfinite(s.B[A0]) && std::isfinite(s.B[D0]));
  std::vector<double> at(na, 0.0), tot(na, 0.0);
  for (size_t i = 0; i < s.segs.size(); ++i) {
    const int a = arm[i];
    if (a < 0 || !std::isfinite(v[i]) || !(w[i] > 0) || !std::isfinite(lev[a])) continue;
    tot[a] += w[i];
    bool base;
    if (cu) base = std::fabs(v[i] - lev[a]) < 0.5;
    else {
      const double x = s.segs[i].mean - lev[a];
      base = have_det && score_component(A0, x, s.E, s.B, s.a1_logistic) < 0.5 && score_component(D0, x, s.E, s.B, s.a1_logistic) < 0.5;
    }
    if (base) at[a] += w[i];
  }
  std::vector<ArmLevel> out(na);
  for (int a = 0; a < na; ++a) out[a] = {lev[a], tot[a] > 0 && (cu || have_det) ? at[a] / tot[a] : NaN};
  return out;
}

// Focal scores: each segment scored on its deviation from its arm's baseline, so a piece at the baseline scores as
// neutral and arm-level change leaves the scores. Segments must be split at the centromere, `base` from
// arm_baselines() on the same segments.  [INFERRED G-14d]
//  - A0, A1 and D0 are scored at the residual. Copy units (a usable ladder, G-13): c' = kmode + (c - b), evaluated at
//    x(max(c', 0)), with x(c) = log2(lad_m * max(lad_r0 + c * step, 1e-6)) the centred log2 value of c copies and
//    kmode = max(1, round(ploidy)) the neutral state the detectors are placed around. Log2 scale (arrays, and
//    sequencing samples flagged step_fallback): mean - b.
//  - D1 is scored at the absolute mean: a homozygous deletion on a lost arm is one copy below its arm (it reads as
//    D0 on the residual), yet it is still a homozygous deletion.
//  - e2 as in score_segments, at the value each component is evaluated at. The sample's own detectors are used.
//  - Where the arm has no baseline (no usable data) the absolute mean is used.
//  - Length rule: contiguous segments of a chromosome in the same direction form a run (gain: residual A0 >= 0.5;
//    loss: residual D0 >= 0.5 or absolute D1 >= 0.5; gain wins a tie). Runs continue across the centromere: the
//    split is a cut of one segment, not a change of state, so a deviation spanning both arms is one event and is
//    measured once against focal_max_bp. A run spanning more than focal_max_bp (first start to
//    last end) is not focal: all four scores and errors of its segments are 0, and its segments' bp are counted in
//    `excluded_frac`, as a share of the sample's segment-covered span (first to last segment per chromosome, summed,
//    as compute_fga).
inline SegScores focal_scores(const Sample& s, const std::vector<Arm>& arms, const std::vector<ArmLevel>& base,
                              const Params& P, double* excluded_frac) {
  const size_t n = s.segs.size();
  const bool cu = uses_copy_units(s);
  const double kmode = std::max(1.0, std::round(s.ploidy));
  auto x = [&](double c) { return std::log2(s.lad_m * std::max(s.lad_r0 + c * s.step, 1e-6)); };
  SegScores out;
  for (int k = 0; k < 4; ++k) { out.s[k].resize(n); out.e2[k].resize(n); }
  for (size_t i = 0; i < n; ++i) {
    const Segment& g = s.segs[i];
    const int a = arm_index(arms, g);
    const double b = a >= 0 && a < (int)base.size() ? base[a].level : NaN;
    double v = g.mean;   // where A0, A1 and D0 are evaluated
    if (std::isfinite(b)) v = cu ? x(std::max(kmode + sample_value(s, g, true) - b, 0.0)) : g.mean - b;
    const double sd = (!is_nan(g.sd) && g.sd > 0) ? g.sd : s.sd_est;
    const double se = sd / std::sqrt((double)std::max(1, g.nmark));
    for (int k = 0; k < 4; ++k) {
      const double at = k == D1 ? g.mean : v;
      out.s[k][i] = (float)score_component(k, at, s.E, s.B, s.a1_logistic);
      const double t = dscore_component(k, at, s.E, s.B, s.a1_logistic) * se;
      out.e2[k][i] = (float)(t * t);
    }
  }
  std::vector<int> dir(n, 0);
  for (size_t i = 0; i < n; ++i)
    dir[i] = out.s[A0][i] >= 0.5f ? 1 : (out.s[D0][i] >= 0.5f || out.s[D1][i] >= 0.5f ? -1 : 0);
  double excluded = 0;
  for (size_t i = 0; i < n;) {
    size_t e = i + 1;
    while (e < n && dir[i] != 0 && dir[e] == dir[i] && s.segs[e].chrom == s.segs[i].chrom) ++e;
    if (dir[i] != 0 && (double)(s.segs[e - 1].bp_end - s.segs[i].bp_start + 1) > P.focal_max_bp)
      for (size_t r = i; r < e; ++r) {
        excluded += (double)(s.segs[r].bp_end - s.segs[r].bp_start + 1);
        for (int k = 0; k < 4; ++k) { out.s[k][r] = 0; out.e2[k][r] = 0; }
      }
    i = e;
  }
  if (excluded_frac) {
    std::vector<int64_t> lo, hi;
    for (const Segment& g : s.segs) {
      if ((int)lo.size() <= g.chrom) { lo.resize(g.chrom + 1, -1); hi.resize(g.chrom + 1, -1); }
      if (lo[g.chrom] < 0 || g.bp_start < lo[g.chrom]) lo[g.chrom] = g.bp_start;
      if (g.bp_end > hi[g.chrom]) hi[g.chrom] = g.bp_end;
    }
    double span = 0;
    for (size_t c = 0; c < lo.size(); ++c) if (lo[c] >= 0) span += (double)(hi[c] - lo[c] + 1);
    *excluded_frac = span > 0 ? excluded / span : NaN;
  }
  return out;
}

// ---- arm-level cohort test  [INFERRED G-14e] ----------------------------------------------------------------
// One row per eligible arm (at least arm_min_markers grid markers), in arm-layout order.
struct ArmResult {
  int arm;                          // index into the arm layout
  int n_samples;                    // samples with an arm baseline (data on the arm)
  double mean[4];                   // cohort means of the four arm scores A0, A1, D0, D1
  double score_gain, score_loss;    // A' = sqrt(A0^2 + A1^2), D' = sqrt(D0^2 + D1^2)
  double p_gain, q_gain, fwer_gain, p_loss, q_loss, fwer_loss;
  double frac_gain, frac_loss;      // share of the n_samples with A0 (D0) >= call_threshold
  double baseline_coverage;         // mean over the samples of the baseline coverage
};

// The four arm scores of a sample: its own detectors at the arm baseline b, on the centred log2 scale. Copy units
// (a usable ladder, G-13): x(b) = log2(lad_m * max(lad_r0 + b * step, 1e-6)); log2 scale: b itself. All four NaN
// where there is no baseline or a score is not finite.
inline void arm_scores(const Sample& s, double b, double out[4]) {
  const double x = uses_copy_units(s) ? std::log2(s.lad_m * std::max(s.lad_r0 + b * s.step, 1e-6)) : b;
  bool ok = std::isfinite(b) && std::isfinite(x);
  for (int k = 0; k < 4 && ok; ++k) { out[k] = score_component(k, x, s.E, s.B, s.a1_logistic); ok = std::isfinite(out[k]); }
  if (!ok) for (int k = 0; k < 4; ++k) out[k] = NaN;
}

// Cohort test of each eligible arm. Observed: per arm, the mean over samples with data of each arm score, combined
// into A' and D'. Null (B = arm_n_perm rounds, one stream Rng(seed, 0x300000000)): each sample's four-score vectors
// are permuted across the eligible arms (Fisher-Yates), keeping its arm burden, and A' and D' are recomputed per
// arm; a sample's missing arm (NaN) moves with the permutation and is skipped, as in permute_block. Gain and loss
// are separate families throughout.
//  - pooled p = (#{null >= obs} + 1) / (N + 1) over all arms and rounds (N = B * n_arm, less any slot with no
//    sample, which gives no value);
//  - q: Benjamini-Hochberg over the arms;
//  - fwer: Westfall-Young max-T, (#{b : max over arms of the null >= obs} + 1) / (B + 1).
// `base` from arm_baselines() per sample (segments split at the centromere); only the detectors and ladder of S
// are used.
inline std::vector<ArmResult> arm_test(const std::vector<Sample>& S, const std::vector<Arm>& arms,
                                       const std::vector<std::vector<ArmLevel>>& base, const Params& P) {
  std::vector<int> el;
  for (int a = 0; a < (int)arms.size(); ++a) if (arms[a].eligible) el.push_back(a);
  const int n = (int)el.size(), N = (int)S.size();
  std::vector<ArmResult> out(n);
  if (n == 0) return out;
  std::vector<double> v((size_t)N * n * 4);   // [sample][eligible arm][component]
  for (int j = 0; j < N; ++j)
    for (int t = 0; t < n; ++t) {
      const int a = el[t];
      arm_scores(S[j], a < (int)base[j].size() ? base[j][a].level : NaN, &v[((size_t)j * n + t) * 4]);
    }
  // A' and D' of every arm for one assignment: arm t takes sample j's vector at slot perm(j, t). cnt 0: no value.
  std::vector<double> acc(4 * (size_t)n);
  std::vector<int> cnt(n);
  auto combine = [&](int t, double& A, double& D) {
    const double c = cnt[t], a0 = acc[4 * t] / c, a1 = acc[4 * t + 1] / c, d0 = acc[4 * t + 2] / c, d1 = acc[4 * t + 3] / c;
    A = std::sqrt(a0 * a0 + a1 * a1); D = std::sqrt(d0 * d0 + d1 * d1);
  };
  auto add = [&](int t, const double* x) {
    if (x[0] != x[0]) return;
    for (int k = 0; k < 4; ++k) acc[4 * t + k] += x[k];
    cnt[t]++;
  };
  // observed
  std::fill(acc.begin(), acc.end(), 0.0); std::fill(cnt.begin(), cnt.end(), 0);
  for (int j = 0; j < N; ++j) for (int t = 0; t < n; ++t) add(t, &v[((size_t)j * n + t) * 4]);
  for (int t = 0; t < n; ++t) {
    ArmResult& r = out[t];
    r.arm = el[t]; r.n_samples = cnt[t];
    for (int k = 0; k < 4; ++k) r.mean[k] = cnt[t] ? acc[4 * t + k] / cnt[t] : NaN;
    r.score_gain = r.score_loss = NaN;
    if (cnt[t]) combine(t, r.score_gain, r.score_loss);
    int fg = 0, fl = 0, nc = 0;
    double cov = 0;
    for (int j = 0; j < N; ++j) {
      const double* x = &v[((size_t)j * n + t) * 4];
      if (x[0] != x[0]) continue;
      fg += x[A0] >= P.call_threshold; fl += x[D0] >= P.call_threshold;
      const double c = el[t] < (int)base[j].size() ? base[j][el[t]].coverage : NaN;
      if (std::isfinite(c)) { cov += c; ++nc; }
    }
    r.frac_gain = cnt[t] ? (double)fg / cnt[t] : NaN;
    r.frac_loss = cnt[t] ? (double)fl / cnt[t] : NaN;
    r.baseline_coverage = nc ? cov / nc : NaN;
  }
  // null
  Rng rng(P.seed, 0x300000000ULL);
  std::vector<int> idx(n);
  for (int t = 0; t < n; ++t) idx[t] = t;
  std::vector<double> nA, nD, tA, tD;
  nA.reserve((size_t)P.arm_n_perm * n); nD.reserve((size_t)P.arm_n_perm * n);
  for (int b = 0; b < P.arm_n_perm; ++b) {
    std::fill(acc.begin(), acc.end(), 0.0); std::fill(cnt.begin(), cnt.end(), 0);
    for (int j = 0; j < N; ++j) {
      for (int t = 0; t + 1 < n; ++t) std::swap(idx[t], idx[t + (int)rng.below((uint32_t)(n - t))]);
      for (int t = 0; t < n; ++t) add(t, &v[((size_t)j * n + idx[t]) * 4]);
    }
    double mA = -1, mD = -1;
    for (int t = 0; t < n; ++t) {
      if (!cnt[t]) continue;
      double A, D;
      combine(t, A, D);
      nA.push_back(A); nD.push_back(D);
      mA = std::max(mA, A); mD = std::max(mD, D);
    }
    if (mA >= 0) { tA.push_back(mA); tD.push_back(mD); }
  }
  for (auto* x : {&nA, &nD, &tA, &tD}) std::sort(x->begin(), x->end());
  auto upper = [](const std::vector<double>& z, double o) {   // (#{z >= o} + 1) / (|z| + 1)
    if (!std::isfinite(o)) return NaN;
    return (double)(z.end() - std::lower_bound(z.begin(), z.end(), o) + 1) / (double)(z.size() + 1);
  };
  std::vector<double> pg, pl;
  std::vector<int> fin;   // arms with an observed score
  for (int t = 0; t < n; ++t) {
    ArmResult& r = out[t];
    r.p_gain = upper(nA, r.score_gain); r.fwer_gain = upper(tA, r.score_gain);
    r.p_loss = upper(nD, r.score_loss); r.fwer_loss = upper(tD, r.score_loss);
    r.q_gain = r.q_loss = NaN;
    if (std::isfinite(r.p_gain)) { fin.push_back(t); pg.push_back(r.p_gain); pl.push_back(r.p_loss); }
  }
  const std::vector<double> qg = bh_adjust(pg), ql = bh_adjust(pl);
  for (size_t i = 0; i < fin.size(); ++i) { out[fin[i]].q_gain = qg[i]; out[fin[i]].q_loss = ql[i]; }
  return out;
}

}  // namespace rae
