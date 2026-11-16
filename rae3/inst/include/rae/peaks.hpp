// RAE 3.0 core -- independence of secondary peaks inside one event (opt-in, peak_independence).  [INFERRED G-14f]
//
// After stage 2 (roi.hpp), an event can hold several peaks. A later peak may only be the shoulder of a higher one:
// the same tumours' alterations run through both. Each later peak is therefore re-scored with the contribution of
// every run (a tumour's contiguous segments at or above a detector threshold) that reaches the apex of an accepted
// peak set to 0, and kept as independent when the residual still reaches the significance cut. The idea is the
// peel-off of GISTIC2 (Mermel et al. 2011 Genome Biol 12:R41), reimplemented from the published description only;
// no GISTIC2 code is used.
//
// Included by pipeline.hpp after Result is defined; include pipeline.hpp, not this file.
#pragma once
#ifndef RAE_RESULT_DEFINED
#error "peaks.hpp needs rae::Result: include pipeline.hpp instead"
#endif
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>
#include "null.hpp"
#include "roi.hpp"
#include "types.hpp"
#include "ubp.hpp"

namespace rae {

// Tumour j's means of components k0 and k0 + 1 over UBP regions r0 .. r0 + nr - 1, computed as score_sample_regions
// computes them but with the scores of the segments flagged in masked[c] set to 0; an empty masked[c] masks nothing,
// and then out equals R.M (test hook rae_masked_means_cpp). out[c][r] is NaN where the tumour has no data.
// [INFERRED G-14f]
inline void masked_region_means(const Result& R, size_t j, int k0, int r0, int nr, const std::vector<char> masked[2],
                                std::vector<float> out[2]) {
  const Profile& pr = R.prof[j];
  const bool wcode = R.P.region_weight_code;
  for (int c = 0; c < 2; ++c) out[c].assign(nr, std::numeric_limits<float>::quiet_NaN());
  size_t pi = 0;
  for (int r = 0; r < nr; ++r) {
    const int32_t rs = R.ubp.start[r0 + r], re = R.ubp.end[r0 + r];
    while (pi < pr.size() && pr[pi].end < rs) ++pi;
    double acc[2] = {0, 0}, w = 0;
    for (size_t z = pi; z < pr.size() && pr[z].start <= re; ++z) {
      const int ov = std::min(pr[z].end, re) - std::max(pr[z].start, rs) + 1;
      if (ov <= 0) continue;
      const double wt = wcode ? (double)ov * ov : (double)ov;   // [CODE D-14]
      for (int c = 0; c < 2; ++c)
        if (masked[c].empty() || !masked[c][pr[z].parent]) acc[c] += wt * (double)R.seg_scores[j].s[k0 + c][pr[z].parent];
      w += wt;
    }
    if (w > 0) for (int c = 0; c < 2; ++c) out[c][r] = (float)(acc[c] / w);
  }
}

// Fills resid_score, resid_p and independent of the peaks of one direction (+1 gain, -1 loss). Rows of the other
// direction, stage-1 rows, primary peaks and peaks alone in their event keep NaN, NaN, -1 (NA).  [INFERRED G-14f]
inline void test_peak_independence(std::vector<Roi>& roi, int direction, const Result& R) {
  const Params& P = R.P;
  const int k0 = direction > 0 ? A0 : D0;           // components: A0, A1 for gain; D0, D1 for loss
  const std::vector<double>& score = direction > 0 ? R.st.Ap : R.st.Dp;
  const std::vector<double>& q = direction > 0 ? R.qA : R.qD;
  const std::vector<double>& nul = direction > 0 ? R.null.A : R.null.D;
  // the cut: the lowest observed score among regions with q <= peak_q (the stage-2 eligibility), genome-wide
  double cut = std::numeric_limits<double>::infinity();
  for (size_t r = 0; r < score.size(); ++r) if (q[r] <= P.peak_q) cut = std::min(cut, score[r]);

  // Runs per tumour and component: maximal stretches of consecutive segments on one chromosome whose score is at
  // least the threshold (peak_run_min; amp_threshold for the linear A1 of the array presets). Under focal_only the
  // segments are the centromere-split ones with focal scores, and a run continues across the centromere: the two
  // halves of one segment, and qualifying neighbours on either arm, are one run.
  struct Run { int32_t s0, s1, m0, m1; };          // segment indices and marker span, inclusive
  const size_t N = R.S.size();
  std::vector<std::vector<Run>> runs[2];
  for (int c = 0; c < 2; ++c) runs[c].resize(N);
  for (size_t j = 0; j < N; ++j) {
    const Sample& s = R.S[j];
    for (int c = 0; c < 2; ++c) {
      const int k = k0 + c;
      const double thr = (k == A1 && !s.a1_logistic) ? P.amp_threshold : P.peak_run_min;
      std::vector<Run>& out = runs[c][j];
      bool open = false;
      for (size_t i = 0; i < s.segs.size(); ++i) {
        const Segment& g = s.segs[i];
        if (!((double)R.seg_scores[j].s[k][i] >= thr)) { open = false; continue; }
        if (open && R.S[j].segs[out.back().s1].chrom == g.chrom) { out.back().s1 = (int32_t)i; out.back().m1 = g.end; }
        else out.push_back({(int32_t)i, (int32_t)i, g.start, g.end});
        open = true;
      }
    }
  }

  // the stage-2 peaks of each event, in order of score (stable: a tie keeps genomic order)
  std::map<int, std::vector<int>> by_event;
  for (size_t x = 0; x < roi.size(); ++x)
    if (roi[x].direction == direction && roi[x].is_peak) by_event[roi[x].event].push_back((int)x);
  for (auto& ev : by_event) {
    std::vector<int>& idx = ev.second;
    if (idx.size() < 2) continue;
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return roi[a].score > roi[b].score; });
    std::vector<int> accepted(1, idx[0]);           // the primary peak: independent stays NA
    for (size_t t = 1; t < idx.size(); ++t) {
      Roi& S = roi[idx[t]];
      const int r0 = S.r_first, nr = S.r_last - S.r_first + 1;
      std::vector<double> sum[2] = {std::vector<double>(nr, 0.0), std::vector<double>(nr, 0.0)};
      std::vector<char> masked[2];
      std::vector<float> mean[2];
      for (size_t j = 0; j < N; ++j) {
        // mask the runs whose marker span meets the apex region of an accepted peak
        bool any = false;
        for (int c = 0; c < 2; ++c) {
          masked[c].assign(R.S[j].segs.size(), 0);
          for (const Run& u : runs[c][j])
            for (int a : accepted) {
              const int lo = R.ubp.start[roi[a].apex], hi = R.ubp.end[roi[a].apex];
              if (u.m0 <= hi && u.m1 >= lo) {
                for (int32_t i = u.s0; i <= u.s1; ++i) masked[c][i] = 1;
                any = true;
                break;
              }
            }
        }
        const size_t base = (size_t)j * R.M.R;
        if (!any) {   // nothing masked: the observed region means, as region_stats adds them
          for (int r = 0; r < nr; ++r) {
            // score_sample_regions marks no data with NaN in all four components at once (w == 0), and region_stats
            // tests only component 0, so testing v[0] is right for the loss components D0, D1 too
            if (R.M.v[0][base + r0 + r] != R.M.v[0][base + r0 + r]) continue;
            for (int c = 0; c < 2; ++c) sum[c][r] += R.M.v[k0 + c][base + r0 + r];
          }
          continue;
        }
        masked_region_means(R, j, k0, r0, nr, masked, mean);
        for (int r = 0; r < nr; ++r) {
          if (mean[0][r] != mean[0][r]) continue;        // no data: not in the cohort mean (as region_stats)
          for (int c = 0; c < 2; ++c) sum[c][r] += mean[c][r];
        }
      }
      // residual summary score per region, over the cohort as before (N unchanged), and its maximum
      double best = 0;
      for (int r = 0; r < nr; ++r) {
        const double n = std::max(1, R.st.n[r0 + r]);
        const double m0 = sum[0][r] / n, m1 = sum[1][r] / n;
        best = std::max(best, std::sqrt(m0 * m0 + m1 * m1));
      }
      S.resid_score = best;
      S.resid_p = empirical_p(std::vector<double>(1, best), nul, P.p_strict)[0];
      S.independent = best >= cut ? 1 : 0;
      if (S.independent) accepted.push_back(idx[t]);
    }
  }
}

}  // namespace rae
