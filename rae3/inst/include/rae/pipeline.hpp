// RAE 3.0 core -- the whole analysis, from parsed input to regions of interest.
#pragma once
#include <map>
#include <string>
#include <vector>
#include "arms.hpp"
#include "density.hpp"
#include "io.hpp"
#include "model.hpp"
#include "null.hpp"
#include "params.hpp"
#include "roi.hpp"
#include "types.hpp"
#include "ubp.hpp"

namespace rae {

struct Result {
  Params P;
  Genome G;
  InputLog log;
  std::vector<Sample> S;
  std::vector<SegScores> seg_scores;
  std::vector<SegScores> abs_scores;  // focal_only: default scores of the split segments, for rae_fga (else empty)
  std::vector<char> member;          // in UBP derivation
  Regions ubp, cleaved_example;      // cleaved_example: first null round, for the size histogram
  ScoreMatrix M;
  RegionStats st;
  NullDist null;
  NullInfo null_info;
  std::vector<double> pA, qA, pD, qD;
  std::vector<Roi> roi;
  double median_segment_bp = NaN;
  // germline (optional normals)  [STATED: M&M Germline]
  int n_normals = 0;
  std::vector<int> germ_gain, germ_loss;  // per UBP region: normals with A0 / D0 >= threshold
  // arm geometry and per-sample arm baselines (focal_only or arm_test; else empty)  [INFERRED G-14c, G-14d]
  std::vector<Arm> arms;
  std::vector<std::vector<ArmLevel>> baselines;   // [sample][arm]
  std::vector<ArmResult> arm_results;             // arm_test: one row per eligible arm  [INFERRED G-14e]
  std::vector<Profile> prof;                      // the observed profiles (one per sample), as scored into M
};
}  // namespace rae
#define RAE_RESULT_DEFINED
#include "peaks.hpp"
namespace rae {

// Steps 1-4 for a set of samples.
inline std::vector<SegScores> fit_samples(std::vector<Sample>& S, const Params& P, uint64_t stream_base) {
  for (Sample& s : S) resolve_sd(s, P, NaN);
  std::vector<double> known;
  for (Sample& s : S) if (!is_nan(s.sd_est)) known.push_back(s.sd_est);
  double cohort = known.empty() ? 0.1 : median(known);
  for (Sample& s : S) if (is_nan(s.sd_est)) s.sd_est = cohort;
  std::vector<SegScores> sc(S.size());
  for (size_t j = 0; j < S.size(); ++j) {
    parameterize_sample(S[j], P, stream_base + j);
    sc[j] = score_segments(S[j]);
  }
  return sc;
}

// Taylor 2010 five-state call (NOT RAE 2008). Precedence [INFERRED G-12c]:
// -2 homozygous deletion, -1 heterozygous loss, 2 amplification, 1 gain, 0 neutral.
inline int discrete_state(double a0, double a1, double d0, double d1, const Params& P) {
  if (d1 >= P.call_threshold) return -2;
  if (d0 >= P.call_threshold) return -1;
  if (a0 >= P.call_threshold) return a1 >= P.amp_threshold ? 2 : 1;
  return 0;
}

// Fraction of the genome altered per sample (NOT RAE 2008), two ways over the same denominator:
//  - rae_fga: segments the RAE scores call altered, A0 or D0 >= burden_threshold (0.75); the CNA burden of
//    Hieronymus 2014. n_altered counts these segments.                       [INFERRED G-12b, G-12e]
//  - fga: segments with log2 ratio (as supplied) > fga_threshold or < -fga_threshold_loss (default: the
//    same 0.2), the cBioPortal FGA; independent of the preset, whose mode-centring can differ.  [INFERRED G-12e]
// Denominator: the sample's segment-covered span (first segment start to last segment end on each
// autosome, summed), as published; or the whole autosomal assembly (fga_assembly).  [INFERRED G-12a]
inline void compute_fga(Result& R) {
  double assembly = 0;
  for (int c = 0; c < R.G.n_chrom; ++c) assembly += (double)R.G.chrom_len[c];
  const double gain_thr = R.P.fga_threshold;
  const double loss_thr = is_nan(R.P.fga_threshold_loss) ? R.P.fga_threshold : std::fabs(R.P.fga_threshold_loss);
  // focal_only: rae_fga and n_altered stay absolute (the default scores of the split segments)  [G-14d]
  const std::vector<SegScores>& sc = R.abs_scores.empty() ? R.seg_scores : R.abs_scores;
  for (size_t j = 0; j < R.S.size(); ++j) {
    Sample& s = R.S[j];
    double rae_bp = 0, thr_bp = 0;
    int n = 0;
    int32_t last_orig = -1;   // focal_only: a centromere-cut pair is one input segment, counted once  [G-14d]
    std::vector<int64_t> lo(R.G.n_chrom, -1), hi(R.G.n_chrom, -1);
    for (size_t i = 0; i < s.segs.size(); ++i) {
      const Segment& g = s.segs[i];
      double len = (double)(g.bp_end - g.bp_start + 1);
      if (sc[j].s[A0][i] >= R.P.burden_threshold || sc[j].s[D0][i] >= R.P.burden_threshold) {
        rae_bp += len;
        if (g.orig < 0 || g.orig != last_orig) ++n;
        last_orig = g.orig;
      }
      const double r = g.mean + s.shift;   // the ratio as supplied
      if (r > gain_thr || r < -loss_thr) thr_bp += len;
      if (lo[g.chrom] < 0 || g.bp_start < lo[g.chrom]) lo[g.chrom] = g.bp_start;
      if (g.bp_end > hi[g.chrom]) hi[g.chrom] = g.bp_end;
    }
    double span = 0;
    for (int c = 0; c < R.G.n_chrom; ++c) if (lo[c] >= 0) span += (double)(hi[c] - lo[c] + 1);
    double denom = R.P.fga_assembly ? assembly : span;
    s.rae_fga = denom > 0 ? rae_bp / denom : NaN;
    s.fga = denom > 0 ? thr_bp / denom : NaN;
    s.n_altered = n;
  }
}

// rows: tumour segments. normal_rows may be empty.
inline Result run(const std::vector<SegRow>& rows, const std::vector<BinRow>& bins,
                  const std::map<std::string, StatRow>& stats, const std::string& build,
                  const Params& P, InputLog log = InputLog(),
                  const std::vector<SegRow>& normal_rows = std::vector<SegRow>(),
                  const std::map<std::string, StatRow>& normal_stats = std::map<std::string, StatRow>(),
                  const std::string& cnv_spec = "default") {
  Result R;
  R.P = P;
  if (P.ploidy_prior && !(P.ploidy_prior_sd > 0 && P.ploidy_prior_wgd >= 0 && P.ploidy_prior_wgd <= 1 && P.ploidy_prior_ambig >= 1))
    throw std::runtime_error("ploidy prior needs ploidy_prior_sd > 0, 0 <= ploidy_prior_wgd <= 1 and ploidy_prior_ambig >= 1");
  if (P.ploidy_prior && !(std::isfinite(P.ploidy_prior_mu1) && std::isfinite(P.ploidy_prior_mu2)))
    throw std::runtime_error("ploidy prior needs finite ploidy_prior_mu1 and ploidy_prior_mu2");
  if (P.focal_only && !(P.focal_max_bp > 0))
    throw std::runtime_error("focal_only needs focal_max_bp > 0");
  if (P.arm_test && !(P.arm_n_perm >= 1 && P.arm_min_markers >= 1))
    throw std::runtime_error("arm_test needs arm_n_perm >= 1 and arm_min_markers >= 1");
  if (P.peak_independence && !(P.peak_run_min > 0 && P.peak_run_min <= 1))
    throw std::runtime_error("peak_independence needs 0 < peak_run_min <= 1");
  if (rows.empty()) throw std::runtime_error("no autosomal segments in input");
  std::vector<SegRow> all(rows);
  all.insert(all.end(), normal_rows.begin(), normal_rows.end());
  R.G = make_genome(build, all, bins, P, log);
  load_cnv(R.G, cnv_spec, log);
  R.S = make_samples(rows, R.G, stats, log);
  if (R.S.empty()) throw std::runtime_error("no usable samples");
  R.seg_scores = fit_samples(R.S, P, 0);                       // steps 1-4
  if (P.focal_only) {
    // [INFERRED G-14d] Detectors, ladder and ploidy come from the unsplit segments (above). Segments are then cut at
    // the centromere, re-scored (the cut changes nmark, hence e2), given arm baselines, and scored on the deviation
    // from them. Everything below (UBP, observed statistic, null, regions, roi, calls, segments) sees the split segments
    // and focal scores; rae_fga, fga and n_altered stay absolute (abs_scores).
    R.arms = arm_layout(R.G, P);
    split_segments_at_centromeres(R.S, R.G, R.arms);
    R.baselines.resize(R.S.size());
    for (size_t j = 0; j < R.S.size(); ++j) {
      R.seg_scores[j] = score_segments(R.S[j]);
      R.baselines[j] = arm_baselines(R.S[j], R.arms, P);
      R.abs_scores.push_back(R.seg_scores[j]);
      R.seg_scores[j] = focal_scores(R.S[j], R.arms, R.baselines[j], P, &R.S[j].focal_excluded_frac);
    }
  }
  R.member = ubp_membership(R.S, P);                           // step 5
  for (size_t j = 0; j < R.S.size(); ++j) if (!R.member[j]) R.S[j].flags |= F_UBP_EXCLUDED;
  std::vector<Profile>& prof = R.prof;
  prof.resize(R.S.size());
  for (size_t j = 0; j < R.S.size(); ++j) prof[j] = profile_of(R.S[j], P);
  R.ubp = build_ubp(prof, R.member, R.G, P);
  R.M = score_matrix(prof, R.seg_scores, R.ubp, true, P.region_weight_code);         // step 6
  R.st = region_stats(R.M, P);                                 // steps 6-7
  R.median_segment_bp = median_segment_bp(R.S, R.G);
  R.null = null_model(R.S, R.seg_scores, R.member, R.G, P, &R.null_info, &R.cleaved_example);  // step 8
  R.pA = empirical_p(R.st.Ap, R.null.A, P.p_strict);                       // step 9, separately for gain and loss
  R.pD = empirical_p(R.st.Dp, R.null.D, P.p_strict);
  R.qA = bh_adjust(R.pA);                                      // [STATED: Main] BH over UBP regions [G-11d]
  R.qD = bh_adjust(R.pD);
  std::vector<Roi> g = find_roi(R.st.Ap, R.st.errA, R.qA, R.ubp.chrom, +1, P);  // step 10
  std::vector<Roi> l = find_roi(R.st.Dp, R.st.errD, R.qD, R.ubp.chrom, -1, P);
  if (P.peak_independence) {                                   // [INFERRED G-14f] opt-in, changes only these rows
    test_peak_independence(g, +1, R);
    test_peak_independence(l, -1, R);
  }
  R.roi = g;
  R.roi.insert(R.roi.end(), l.begin(), l.end());
  compute_fga(R);
  // step 11: germline from normals, on UBP regions  [INFERRED G-11g, G-11h]
  R.germ_gain.assign(R.ubp.size(), 0); R.germ_loss.assign(R.ubp.size(), 0);
  if (!normal_rows.empty()) {
    InputLog nlog;
    std::vector<Sample> Nn = make_samples(normal_rows, R.G, normal_stats, nlog);
    std::vector<SegScores> nsc = fit_samples(Nn, P, 0x200000000ULL);
    std::vector<Profile> np(Nn.size());
    for (size_t j = 0; j < Nn.size(); ++j) np[j] = profile_of(Nn[j], P);
    ScoreMatrix NM = score_matrix(np, nsc, R.ubp, false);
    R.n_normals = (int)Nn.size();
    for (int j = 0; j < NM.N; ++j)
      for (int r = 0; r < NM.R; ++r) {
        if (NM.at(A0, r, j) >= P.germline_min_score) R.germ_gain[r]++;
        if (NM.at(D0, r, j) >= P.germline_min_score) R.germ_loss[r]++;
      }
  }
  // Arm-level cohort test, last; it changes nothing above.  [INFERRED G-14e] With focal_only the arms and baselines
  // are those of the focal analysis; otherwise they come from a copy of the samples split at the centromere, so the
  // main analysis keeps its unsplit segments.
  if (P.arm_test) {
    if (!P.focal_only) {
      R.arms = arm_layout(R.G, P);
      std::vector<Sample> T(R.S);
      split_segments_at_centromeres(T, R.G, R.arms);
      R.baselines.resize(T.size());
      for (size_t j = 0; j < T.size(); ++j) R.baselines[j] = arm_baselines(T[j], R.arms, P);
    }
    R.arm_results = arm_test(R.S, R.arms, R.baselines, P);
  }
  // one warning for both switches when the build has no centromere coordinates  [INFERRED G-14d, G-14e]
  for (const Arm& a : R.arms)
    if (a.arm == 'w') {
      const std::string who = P.focal_only && P.arm_test ? "focal_only, arm_test" : P.focal_only ? "focal_only" : "arm_test";
      log.warnings.push_back(who + ": no centromere coordinates for build '" + R.G.build + "'; each chromosome is one arm");
      break;
    }
  R.log = log;
  return R;
}

}  // namespace rae
