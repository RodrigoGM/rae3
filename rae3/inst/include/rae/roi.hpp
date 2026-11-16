// RAE 3.0 core -- regions of interest, spec step 10.
#pragma once
#include <algorithm>
#include <vector>
#include "density.hpp"
#include "types.hpp"

namespace rae {

struct Roi {
  int direction = 0;      // +1 gain, -1 loss
  int event = 0;          // stage-1 event id (within direction)
  bool is_peak = false;   // false: stage-1 region; true: stage-2 peak inside it
  int r_first = 0, r_last = 0;  // UBP region indices, inclusive
  int apex = -1;          // UBP region of maximum score
  double score = 0, error = 0, q = 1;  // at the apex
  int n_regions = 0;
  // independence of a secondary peak (peak_independence only; peaks.hpp)  [INFERRED G-14f]
  double resid_score = NaN, resid_p = NaN;
  int independent = -1;   // -1 NA (stage-1 row, primary peak, or the only peak of its event), 0 no, 1 yes
};

// score, err, q: per UBP region for one direction. chrom: region chromosome.
// [CODE D-13] Stage-2 peak detection exactly as coded in getStageTwoPeaks / detectCopyNumberPeaks.
//  - the score trace is floored at cut2, the lowest score among regions with q <= peak_q genome-wide,
//    and padded with cut2 at both ends;
//  - errors are divided by roi_error_factor (0.5, i.e. doubled);
//  - a running-extremum detector reports a maximum once the trace has fallen by
//    shoulder x mean(error) below it;
//  - a peak is dropped if score - error < cut2;
//  - neighbours are merged while contiguous and within +/- error of the peak score (two-sided).
inline void coded_peaks(const std::vector<double>& score, const std::vector<double>& err,
                        const std::vector<double>& q, int i, int j, double cut2, int direction, int ev,
                        const Params& P, std::vector<Roi>& out) {
  const int n = j - i + 1;
  std::vector<double> trace(n + 2, cut2), delta(n);
  double mean_delta = 0;
  for (int t = 0; t < n; ++t) {
    trace[t + 1] = std::max(score[i + t], cut2);
    delta[t] = err[i + t] / P.roi_error_factor;
    mean_delta += delta[t];
  }
  mean_delta /= n;
  const double thr = mean_delta * P.shoulder_factor;
  std::vector<int> pks;
  bool lookformax = true;
  double mx = -1e300, mn = 1e300;
  int mxpos = -1;
  for (int t = 0; t < n + 2; ++t) {
    double v = trace[t];
    if (v > mx) { mx = v; mxpos = t; }
    if (v < mn) mn = v;
    if (lookformax) { if (v < mx - thr) { pks.push_back(mxpos); mn = v; lookformax = false; } }
    else if (v > mn + thr) { mx = v; mxpos = t; lookformax = true; }
  }
  std::vector<char> merged_flag(n, 0);
  for (int pos : pks) {
    int p = pos - 1;  // index within the event
    if (p < 0 || p >= n) continue;
    if (score[i + p] - delta[p] < cut2) continue;
    if (merged_flag[p]) continue;
    int a = p, b = p;
    double minq = q[i + p];
    while (b + 1 < n && std::fabs(score[i + b + 1] - score[i + p]) <= delta[p]) { ++b; merged_flag[b] = 1; minq = std::min(minq, q[i + b]); }
    while (a - 1 >= 0 && std::fabs(score[i + a - 1] - score[i + p]) <= delta[p]) { --a; merged_flag[a] = 1; minq = std::min(minq, q[i + a]); }
    Roi r;
    r.direction = direction; r.event = ev; r.is_peak = true; r.r_first = i + a; r.r_last = i + b;
    r.apex = i + p; r.score = score[i + p]; r.error = err[i + p]; r.q = minq; r.n_regions = b - a + 1;
    out.push_back(r);
  }
}

inline std::vector<Roi> find_roi(const std::vector<double>& score, const std::vector<double>& err,
                                 const std::vector<double>& q, const std::vector<int32_t>& chrom,
                                 int direction, const Params& P) {
  std::vector<Roi> out;
  const int R = (int)score.size();
  int ev = 0;
  double cut2 = 1e300;  // [CODE D-13]
  std::vector<int> sig_on_chrom;
  if (P.roi_code) {
    for (int r = 0; r < R; ++r) {
      if (q[r] <= P.peak_q) cut2 = std::min(cut2, score[r]);
      if (q[r] <= P.fdr) {
        if ((int)sig_on_chrom.size() <= chrom[r]) sig_on_chrom.resize(chrom[r] + 1, 0);
        sig_on_chrom[chrom[r]]++;
      }
    }
  }
  for (int i = 0; i < R;) {
    if (!(q[i] <= P.fdr)) { ++i; continue; }
    // Stage 1: isolated significant regions, or physically adjacent ones merged to the largest
    // boundaries of the event.  [STATED: Main ROI]
    int j = i;
    while (j + 1 < R && chrom[j + 1] == chrom[i] && q[j + 1] <= P.fdr) ++j;
    Roi e;
    e.direction = direction; e.event = ev; e.r_first = i; e.r_last = j; e.n_regions = j - i + 1;
    e.apex = i;
    for (int r = i; r <= j; ++r) if (score[r] > score[e.apex]) e.apex = r;
    e.score = score[e.apex]; e.error = err[e.apex]; e.q = q[e.apex];
    out.push_back(e);
    if (P.roi_code) {
      // As coded: peaks are sought in any event containing a region with q <= peak_q, unless it is
      // the only significant region on its chromosome.
      bool any = false;
      for (int r = i; r <= j; ++r) any = any || q[r] <= P.peak_q;
      if (any && sig_on_chrom[chrom[i]] > 1) coded_peaks(score, err, q, i, j, cut2, direction, ev, P, out);
    } else if (j > i) {
      // Stage 2  [STATED: M&M ROI]: "in each of these larger regions spanning more than one region"
      // shoulder sensitivity = 2 x median analytical error of all UBP regions in the event  [G-09f]
      std::vector<double> ev_err(err.begin() + i, err.begin() + j + 1);
      const double shoulder = P.shoulder_factor * median(ev_err);
      std::vector<char> taken(j - i + 1, 0);
      for (int r = i; r <= j; ++r) {
        if (!(q[r] <= P.peak_q)) continue;  // only regions an order of magnitude more significant
        // local maximum of the summary score; a plateau is one peak with its 5'-most region as apex [G-09b]
        if (r > i && score[r - 1] >= score[r]) continue;
        int pe = r;
        while (pe + 1 <= j && score[pe + 1] == score[r]) ++pe;
        if (pe + 1 <= j && score[pe + 1] > score[r]) continue;
        // prominence: drop between this maximum and any higher region, or the event edge  [INFERRED G-09a]
        double minL = score[r], minR = score[r];
        bool higherL = false, higherR = false;
        for (int t = r - 1; t >= i; --t) { if (score[t] > score[r]) { higherL = true; break; } minL = std::min(minL, score[t]); }
        for (int t = pe + 1; t <= j; ++t) { if (score[t] > score[r]) { higherR = true; break; } minR = std::min(minR, score[t]); }
        double ref;
        if (higherL && higherR) ref = std::max(minL, minR);
        else if (higherL) ref = minL;
        else if (higherR) ref = minR;
        else ref = std::min(minL, minR);  // the event's highest point: measured against the event floor
        if (score[r] - ref < shoulder || score[r] - ref <= 0) continue;
        // "Regions adjacent to a peak are merged with that locus if their summary score fell within
        // the error bar of the peak"  [STATED: M&M ROI]; contiguous, both directions, inside the
        // event  [INFERRED G-09c, G-09d]
        int a = r, b = pe;
        while (a - 1 >= i && score[a - 1] >= score[r] - err[r]) --a;
        while (b + 1 <= j && score[b + 1] >= score[r] - err[r]) ++b;
        Roi p;
        p.direction = direction; p.event = ev; p.is_peak = true; p.r_first = a; p.r_last = b;
        p.apex = r; p.score = score[r]; p.error = err[r]; p.q = q[r]; p.n_regions = b - a + 1;
        out.push_back(p);
        r = pe;
      }
    }
    ++ev;
    i = j + 1;
  }
  return out;
}

}  // namespace rae
