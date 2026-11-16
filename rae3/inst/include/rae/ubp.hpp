// RAE 3.0 core -- unified breakpoint profile (spec step 5) and region scoring (steps 6-7).
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "density.hpp"
#include "model.hpp"
#include "types.hpp"

namespace rae {

// A piece of a sample's profile: marker range plus the index of the original segment whose value it
// carries. Original segments are pieces of themselves; cleaved segments inherit their parent's value.
// [STATED: M&M Null] "Each split inherited the copy number assigned to its parental segment"
struct Piece { int32_t start, end, parent; bool ubp_start; };
typedef std::vector<Piece> Profile;

inline Profile profile_of(const Sample& s, const Params& P) {
  Profile p(s.segs.size());
  for (size_t i = 0; i < s.segs.size(); ++i)
    p[i] = {s.segs[i].start, s.segs[i].end, (int32_t)i,
            s.segs[i].nmark >= P.ubp_min_markers};  // [STATED: M&M UBP] "segments of >= 3 markers" [G-07f]
  return p;
}

// Samples whose segment count is an outlier are left out of UBP derivation (but are still scored).
// [STATED: M&M UBP]; the criterion is not published: Tukey fences on log10(count). [INFERRED G-07a]
inline std::vector<char> ubp_membership(const std::vector<Sample>& S, const Params& P) {
  std::vector<char> in(S.size(), 1);
  if (P.ubp_outlier_k <= 0 || S.size() < 10) return in;
  // Segment counts. The original code counts segments of >= 3 markers.  [CODE D-09]
  std::vector<double> lc(S.size());
  for (size_t j = 0; j < S.size(); ++j) {
    double c = 0;
    if (P.ubp_outlier_log) c = (double)S[j].segs.size();
    else for (const Segment& g : S[j].segs) c += g.nmark >= P.ubp_min_markers;
    lc[j] = P.ubp_outlier_log ? std::log10(std::max(1.0, c)) : c;
  }
  std::vector<double> t(lc);
  std::sort(t.begin(), t.end());
  double q1, q3;
  if (P.ubp_outlier_log) { q1 = quantile7_sorted(t, 0.25); q3 = quantile7_sorted(t, 0.75); }
  else {  // R boxplot.stats: Tukey hinges (fivenum), not type-7 quartiles  [CODE D-09]
    size_t n = t.size();
    double n4 = std::floor((n + 3) / 2.0) / 2.0;
    auto at = [&](double pos) { return 0.5 * (t[(size_t)std::floor(pos) - 1] + t[(size_t)std::ceil(pos) - 1]); };
    q1 = at(n4); q3 = at(n + 1 - n4);
  }
  const double iqr = q3 - q1;
  for (size_t j = 0; j < S.size(); ++j)
    if (lc[j] < q1 - P.ubp_outlier_k * iqr || lc[j] > q3 + P.ubp_outlier_k * iqr) in[j] = 0;
  return in;
}

// Build a breakpoint profile from the unique segment starts of the member samples.
inline Regions build_ubp(const std::vector<Profile>& prof, const std::vector<char>& member,
                         const Genome& G, const Params& P) {
  const int M = G.n_markers();
  std::vector<char> is_start(M, 0);
  for (int c = 0; c < G.n_chrom; ++c)
    if (G.chrom_first[c] < G.chrom_first[c + 1]) is_start[G.chrom_first[c]] = 1;  // [INFERRED G-07c]
  if (P.split_at_centromere)  // [INFERRED G-07e] off by default
    for (int c = 0; c < G.n_chrom; ++c)
      for (int m = G.chrom_first[c]; m < G.chrom_first[c + 1]; ++m)
        if (G.cen_end[c] > 0 && G.m_start[m] > G.cen_end[c]) { is_start[m] = 1; break; }
  for (size_t j = 0; j < prof.size(); ++j) {
    if (!member[j]) continue;
    for (const Piece& p : prof[j]) if (p.ubp_start) is_start[p.start] = 1;  // [STATED: M&M UBP] unique starts
  }
  // Regions: each start runs to the marker directly 5' of the next start.  [STATED: M&M UBP]
  Regions raw;
  for (int c = 0; c < G.n_chrom; ++c) {
    int lo = G.chrom_first[c], hi = G.chrom_first[c + 1];
    int cur = -1;
    for (int m = lo; m < hi; ++m) {
      if (is_start[m]) {
        if (cur >= 0) { raw.start.push_back(cur); raw.end.push_back(m - 1); raw.chrom.push_back(c); }
        cur = m;
      }
    }
    if (cur >= 0) { raw.start.push_back(cur); raw.end.push_back(hi - 1); raw.chrom.push_back(c); }
  }
  // Single-marker regions.  [STATED: M&M UBP]
  //  (ii) a contiguous stretch of single-marker regions is consolidated as one independent region;
  //  (i)  an isolated single-marker region is merged into the physically adjacent 3' region
  //       (into the 5' region when it is the last one on its chromosome  [INFERRED G-07b]).
  Regions cons;
  const int R = raw.size();
  for (int i = 0; i < R;) {
    bool single = raw.start[i] == raw.end[i];
    int j = i;
    if (single) while (j + 1 < R && raw.chrom[j + 1] == raw.chrom[i] && raw.start[j + 1] == raw.end[j + 1]) ++j;
    cons.start.push_back(raw.start[i]); cons.end.push_back(raw.end[j]); cons.chrom.push_back(raw.chrom[i]);
    i = j + 1;
  }
  Regions out;
  const int C = cons.size();
  for (int i = 0; i < C; ++i) {
    bool single = cons.start[i] == cons.end[i];
    if (single && i + 1 < C && cons.chrom[i + 1] == cons.chrom[i]) {
      cons.start[i + 1] = cons.start[i];  // merge into the 3' neighbour
      continue;
    }
    if (single && !out.start.empty() && out.chrom.back() == cons.chrom[i]) {
      out.end.back() = cons.end[i];       // last on chromosome: merge 5'
      continue;
    }
    out.start.push_back(cons.start[i]); out.end.push_back(cons.end[i]); out.chrom.push_back(cons.chrom[i]);
  }
  return out;
}

// Region x sample matrix of component scores, sample-major: v[k][j * R + r].
// Missing (sample has no data on that chromosome) is NaN.
struct ScoreMatrix {
  int R = 0, N = 0;
  std::vector<float> v[4];
  std::vector<float> e2[4];  // squared error terms; empty when not requested
  float at(int k, int r, int j) const { return v[k][(size_t)j * R + r]; }
};

// [STATED: M&M] "the weighted arithmetic mean for each S_k across the m markers in that region";
// weights are the marker counts of each overlapping piece.  [INFERRED G-06a, G-10e]
inline void score_sample_regions(const Profile& prof, const SegScores& sc, const Regions& reg,
                                 ScoreMatrix& M, int j, bool with_error, bool weight_code = false) {
  const int R = reg.size();
  size_t base = (size_t)j * R;
  size_t pi = 0;
  const size_t np = prof.size();
  for (int r = 0; r < R; ++r) {
    while (pi < np && prof[pi].end < reg.start[r]) ++pi;
    double acc[4] = {0, 0, 0, 0}, acc_e[4] = {0, 0, 0, 0}, w = 0;
    for (size_t q = pi; q < np && prof[q].start <= reg.end[r]; ++q) {
      int ov = std::min(prof[q].end, reg.end[r]) - std::max(prof[q].start, reg.start[r]) + 1;
      if (ov <= 0) continue;
      // [CODE D-14] the code gives each marker the weight (run length / region length), so a run of
      // g markers contributes with weight g^2; the paper's marker mean is weight g.
      const double wt = weight_code ? (double)ov * ov : (double)ov;
      for (int k = 0; k < 4; ++k) {
        acc[k] += wt * (double)sc.s[k][prof[q].parent];
        if (with_error) acc_e[k] += wt * (double)sc.e2[k][prof[q].parent];
      }
      w += wt;
    }
    for (int k = 0; k < 4; ++k) {
      M.v[k][base + r] = w > 0 ? (float)(acc[k] / w) : std::numeric_limits<float>::quiet_NaN();
      if (with_error) M.e2[k][base + r] = w > 0 ? (float)(acc_e[k] / w) : 0.f;
    }
  }
}

inline ScoreMatrix score_matrix(const std::vector<Profile>& prof, const std::vector<SegScores>& sc,
                                const Regions& reg, bool with_error, bool weight_code = false) {
  ScoreMatrix M;
  M.R = reg.size(); M.N = (int)prof.size();
  for (int k = 0; k < 4; ++k) {
    M.v[k].resize((size_t)M.R * M.N);
    if (with_error) M.e2[k].resize((size_t)M.R * M.N);
  }
  for (int j = 0; j < M.N; ++j) score_sample_regions(prof[j], sc[j], reg, M, j, with_error, weight_code);
  return M;
}

// Cohort-level region statistics: component means, summary scores (Eq. 3) and analytical error.
struct RegionStats {
  std::vector<double> mean[4], err[4];  // cohort mean of S_k and delta_k
  std::vector<double> Ap, Dp, errA, errD;
  std::vector<int> n;                   // samples with data
};

inline RegionStats region_stats(const ScoreMatrix& M, const Params& P) {
  RegionStats st;
  for (int k = 0; k < 4; ++k) { st.mean[k].assign(M.R, 0); st.err[k].assign(M.R, 0); }
  st.n.assign(M.R, 0);
  for (int j = 0; j < M.N; ++j) {
    size_t base = (size_t)j * M.R;
    for (int r = 0; r < M.R; ++r) {
      float v0 = M.v[0][base + r];
      if (v0 != v0) continue;
      st.n[r]++;
      for (int k = 0; k < 4; ++k) {
        st.mean[k][r] += M.v[k][base + r];
        if (!M.e2[k].empty()) st.err[k][r] += M.e2[k][base + r];
      }
    }
  }
  st.Ap.resize(M.R); st.Dp.resize(M.R); st.errA.resize(M.R); st.errD.resize(M.R);
  for (int r = 0; r < M.R; ++r) {
    double N = std::max(1, st.n[r]);
    for (int k = 0; k < 4; ++k) {
      st.mean[k][r] /= N;                      // [STATED: M&M] "the mean S_k across all samples"
      // [STATED: S1-C as rendered] delta_k = sqrt(sum_j (dS s/sqrt(n))^2 / N)   [INFERRED G-10b]
      st.err[k][r] = P.error_sem ? std::sqrt(st.err[k][r]) / N : std::sqrt(st.err[k][r] / N);
    }
    // [STATED: Eq. 3]  A' = sqrt(A0^2 + A1^2),  D' = sqrt(D0^2 + D1^2)
    double a0 = st.mean[A0][r], a1 = st.mean[A1][r], d0 = st.mean[D0][r], d1 = st.mean[D1][r];
    st.Ap[r] = std::sqrt(a0 * a0 + a1 * a1);
    st.Dp[r] = std::sqrt(d0 * d0 + d1 * d1);
    // [STATED: S1-C]  delta_A = sqrt((A0^2 dA0^2 + A1^2 dA1^2) / A'^2)
    auto comb = [](double x0, double e0, double x1, double e1, double norm) {
      return norm > 0 ? std::sqrt((x0 * x0 * e0 * e0 + x1 * x1 * e1 * e1) / (norm * norm)) : 0.0;
    };
    st.errA[r] = comb(a0, st.err[A0][r], a1, st.err[A1][r], st.Ap[r]);
    st.errD[r] = comb(d0, st.err[D0][r], d1, st.err[D1][r], st.Dp[r]);
  }
  return st;
}

}  // namespace rae
