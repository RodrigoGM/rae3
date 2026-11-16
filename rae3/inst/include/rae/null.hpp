// RAE 3.0 core -- background (null) model, spec step 8, and p-values, spec step 9.
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "density.hpp"
#include "rng.hpp"
#include "types.hpp"
#include "ubp.hpp"

namespace rae {

// Randomised recombination hotspots, as marker indices at which a cut may be placed.
// [STATED: M&M Null] per arm: R(h) = hotspot distribution "with a width equal to the median distance
// between all hotspots on that arm", "normalized to range between 0 and 1"; draw x ~ U(arm),
// u ~ U[0,1]; accept if u < R(x); until the accepted count equals the real count on that arm.
// R(h) is a histogram with that bin width.  [INFERRED G-08a]
// Without a hotspot table: uniform positions at a fixed mean spacing, a deviation.  [INFERRED G-08k]
inline std::vector<int32_t> randomized_cuts(const Genome& G, const Params& P, Rng& rng) {
  std::vector<int32_t> cuts;
  for (int c = 0; c < G.n_chrom; ++c) {
    int lo = G.chrom_first[c], hi = G.chrom_first[c + 1];
    if (hi - lo < 2) continue;
    const int64_t clen = G.chrom_len[c] > 0 ? G.chrom_len[c] : G.m_end[hi - 1];
    std::vector<int64_t> pos;
    for (int arm = 0; arm < 2; ++arm) {
      int64_t a0 = arm == 0 ? 1 : (G.cen_end[c] > 0 ? G.cen_end[c] + 1 : 1);
      int64_t a1 = arm == 0 ? (G.cen_start[c] > 0 ? G.cen_start[c] : clen) : clen;
      if (arm == 1 && G.cen_end[c] <= 0) break;  // unknown centromere: whole chromosome is one "arm"
      if (a1 <= a0) continue;
      if (G.has_hotspots) {
        std::vector<int64_t> h;
        for (int64_t x : G.hotspots[c]) if (x >= a0 && x <= a1) h.push_back(x);
        if (h.size() < 2) continue;  // [STATED: M&M Null] non-acrocentric arms only  [G-08b]
        std::vector<double> gaps;
        for (size_t i = 1; i < h.size(); ++i) gaps.push_back((double)(h[i] - h[i - 1]));
        double w = std::max(1.0, median(gaps));
        std::vector<int64_t> hits;
        size_t accepted = 0;
        if (P.hotspot_kde) {
          // [CODE D-10] R density(hot, bw = width): Gaussian KDE on 512 points over
          // [min - 3 bw, max + 3 bw], scaled to a maximum of 1; zero probability outside that range.
          const int ng = 512;
          const double lo_x = (double)h.front() - 3 * w, hi_x = (double)h.back() + 3 * w, dx = (hi_x - lo_x) / (ng - 1);
          std::vector<double> kd(ng, 0.0);
          for (int64_t x : h) {
            int c0 = (int)(((double)x - lo_x) / dx);
            int span = (int)std::ceil(5 * w / dx) + 1;
            for (int g = std::max(0, c0 - span); g <= std::min(ng - 1, c0 + span); ++g) {
              double z = (lo_x + g * dx - (double)x) / w;
              kd[g] += std::exp(-0.5 * z * z);
            }
          }
          double mx = *std::max_element(kd.begin(), kd.end());
          while (accepted < h.size()) {
            double x = a0 + rng.unif() * (double)(a1 - a0);
            double u = rng.unif();
            double pos = (x - lo_x) / dx;
            double lev = -1;
            if (pos >= 0 && pos <= ng - 1) {
              int g = std::min(ng - 2, (int)pos);
              lev = (kd[g] + (pos - g) * (kd[g + 1] - kd[g])) / mx;
            }
            if (u <= lev) { ++accepted; hits.push_back((int64_t)x); }
          }
        } else {
          size_t nb = (size_t)std::floor((a1 - a0) / w) + 1;
          std::vector<double> hist(nb, 0.0);
          for (int64_t x : h) hist[(size_t)((x - a0) / w)] += 1;
          double mx = *std::max_element(hist.begin(), hist.end());
          while (accepted < h.size()) {
            double x = a0 + rng.unif() * (double)(a1 - a0);
            double u = rng.unif();
            if (u < hist[(size_t)((x - a0) / w)] / mx) { hits.push_back((int64_t)x); ++accepted; }
          }
        }
        pos.insert(pos.end(), hits.begin(), hits.end());
      } else {
        size_t n = (size_t)((a1 - a0) / P.uniform_cleave_spacing);
        for (size_t i = 0; i < n; ++i) pos.push_back(a0 + (int64_t)(rng.unif() * (double)(a1 - a0)));
      }
    }
    // "mapped to the marker nearest each in the genome"  [STATED: M&M Null]
    for (int64_t x : pos) {
      int a = lo, b = hi - 1;
      while (a < b) { int m = (a + b) / 2; if ((G.m_start[m] + G.m_end[m]) / 2 < x) a = m + 1; else b = m; }
      int best = a;
      if (a > lo) {
        int64_t da = std::llabs((G.m_start[a] + G.m_end[a]) / 2 - x);
        int64_t db = std::llabs((G.m_start[a - 1] + G.m_end[a - 1]) / 2 - x);
        if (db < da) best = a - 1;
      }
      cuts.push_back(best);
    }
  }
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  return cuts;
}

inline int64_t bp_span(const Genome& G, int32_t a, int32_t b) { return G.m_end[b] - G.m_start[a] + 1; }

// Median size (bp) of all segments in the cohort.  [STATED: M&M Null] "median size of all segments by
// disease type"; base pairs  [INFERRED G-08c]
inline double median_segment_bp(const std::vector<Sample>& S, const Genome& G) {
  std::vector<double> v;
  for (const Sample& s : S) for (const Segment& g : s.segs) v.push_back((double)bp_span(G, g.start, g.end));
  return median(v);
}

// Cleave one sample at the randomised cut markers.  [STATED: M&M Null]
//  - only segments larger than the cohort median size are cut;
//  - both resulting pieces must have >= 2 markers and be >= half the median size;
//  - pieces inherit the parent's value.
// Several cuts inside one segment are applied greedily 5' -> 3'.  [INFERRED G-08d]
inline Profile cleave(const Sample& s, const std::vector<int32_t>& cuts, const Genome& G,
                      double med_bp, const Params& P) {
  Profile out;
  out.reserve(s.segs.size() * 2);
  for (size_t i = 0; i < s.segs.size(); ++i) {
    const Segment& g = s.segs[i];
    int32_t last = g.start;
    if ((double)bp_span(G, g.start, g.end) > med_bp) {
      auto it = std::upper_bound(cuts.begin(), cuts.end(), g.start);
      if (P.cleave_code) {
        // [CODE D-11] As coded: walking the hotspot markers inside the segment 5' -> 3', a piece is
        // closed at a hotspot when the distance from the piece's start exceeds half the median size.
        // The last hotspot in a segment never closes a piece (the piece runs to the segment end), and
        // there is no minimum-marker test and no test on the remainder.
        std::vector<int32_t> spots;
        auto lo_it = std::lower_bound(cuts.begin(), cuts.end(), g.start);
        for (auto jt = lo_it; jt != cuts.end() && *jt <= g.end; ++jt) spots.push_back(*jt);
        for (size_t k = 0; k + 1 < spots.size(); ++k) {
          double size = (double)(G.m_start[spots[k]] - G.m_start[last]);
          if (size > 0.5 * med_bp && spots[k] < g.end) {
            out.push_back({last, spots[k], (int32_t)i, false});
            last = spots[k] + 1;
          }
        }
      } else {
        for (; it != cuts.end() && *it <= g.end; ++it) {
          int32_t c = *it;  // a new piece would start at marker c
          if (c - last < 2 || g.end - c + 1 < 2) continue;
          if ((double)bp_span(G, last, c - 1) < 0.5 * med_bp) continue;
          if ((double)bp_span(G, c, g.end) < 0.5 * med_bp) continue;
          out.push_back({last, c - 1, (int32_t)i, false});
          last = c;
        }
      }
    }
    out.push_back({last, g.end, (int32_t)i, false});
  }
  // The >= 3 marker rule for UBP starts applies to the cleaved pieces as it did to segments. A piece
  // that begins where its parent begins keeps the parent's eligibility only if it is still >= 3 markers.
  for (Piece& p : out) {
    int span = p.end - p.start + 1;
    const Segment& par = s.segs[p.parent];
    int est = (p.start == par.start && p.end == par.end) ? par.nmark : span;
    p.ubp_start = est >= P.ubp_min_markers;
  }
  return out;
}

struct NullDist { std::vector<double> A, D; };  // pooled null summary scores, sorted ascending

// One block of permutations on one cleaved breakpoint profile.
// [STATED: Main] "the values of the multi-component model in each region of this UBP (A0, A1, D0, D1)
// are permuted together to another position of the UBP in each sample and recombined across tumors".
// Regions are exchangeable units [INFERRED G-08e]; the first `null_positions` positions of each
// permutation are scored [INFERRED G-08f]; all samples take part [INFERRED G-08j].
inline void permute_block(const ScoreMatrix& M, int n_perm, const Params& P, Rng& rng, NullDist& out) {
  const int R = M.R, N = M.N;
  if (R == 0 || N == 0) return;
  const int K = (P.null_positions <= 0 || P.null_positions > R) ? R : P.null_positions;
  std::vector<int32_t> idx(R);
  for (int r = 0; r < R; ++r) idx[r] = r;
  std::vector<double> acc(4 * (size_t)K);
  std::vector<int> cnt(K);
  for (int b = 0; b < n_perm; ++b) {
    std::fill(acc.begin(), acc.end(), 0.0);
    std::fill(cnt.begin(), cnt.end(), 0);
    for (int j = 0; j < N; ++j) {
      const size_t base = (size_t)j * R;
      for (int t = 0; t < K; ++t) {  // partial Fisher-Yates: idx[0..K) is a random K-subset in random order
        int u = t + (int)rng.below((uint32_t)(R - t));
        std::swap(idx[t], idx[u]);
        int r = idx[t];
        float v0 = M.v[0][base + r];
        if (v0 != v0) continue;
        acc[4 * t + 0] += v0;
        acc[4 * t + 1] += M.v[1][base + r];
        acc[4 * t + 2] += M.v[2][base + r];
        acc[4 * t + 3] += M.v[3][base + r];
        cnt[t]++;
      }
    }
    for (int t = 0; t < K; ++t) {
      double n = std::max(1, cnt[t]);
      double a0 = acc[4 * t] / n, a1 = acc[4 * t + 1] / n, d0 = acc[4 * t + 2] / n, d1 = acc[4 * t + 3] / n;
      out.A.push_back(std::sqrt(a0 * a0 + a1 * a1));
      out.D.push_back(std::sqrt(d0 * d0 + d1 * d1));
    }
  }
}

struct NullInfo { int rounds = 0; std::vector<int> cleaved_regions; std::vector<int> n_cuts;
                  std::vector<int32_t> first_cuts; };  // first_cuts: marker indices of round 0 (diagnostic plots)

// Full null model: hotspot randomisation, cleavage, UBP regeneration and re-scoring are redone
// every `recleave_every` permutations.  [STATED: M&M Null] reading [INFERRED G-08g]
inline NullDist null_model(const std::vector<Sample>& S, const std::vector<SegScores>& sc,
                           const std::vector<char>& member, const Genome& G, const Params& P,
                           NullInfo* info = nullptr, Regions* first_cleaved = nullptr) {
  NullDist nd;
  const double med_bp = median_segment_bp(S, G);
  int done = 0, round = 0;
  while (done < P.n_perm) {
    int block = std::min(P.recleave_every > 0 ? P.recleave_every : P.n_perm, P.n_perm - done);
    Rng rng(P.seed, 0x100000000ULL + (uint64_t)round);
    std::vector<int32_t> cuts = randomized_cuts(G, P, rng);
    std::vector<Profile> prof(S.size());
    for (size_t j = 0; j < S.size(); ++j) prof[j] = cleave(S[j], cuts, G, med_bp, P);
    Regions reg = build_ubp(prof, member, G, P);  // "the UBP is derived again"  [STATED: Main]
    ScoreMatrix M = score_matrix(prof, sc, reg, false, P.region_weight_code);
    permute_block(M, block, P, rng, nd);
    if (info) { info->cleaved_regions.push_back(reg.size()); info->n_cuts.push_back((int)cuts.size()); }
    if (first_cleaved && round == 0) *first_cleaved = reg;
    if (info && round == 0) info->first_cuts = cuts;
    done += block; ++round;
  }
  if (info) info->rounds = round;
  std::sort(nd.A.begin(), nd.A.end());
  std::sort(nd.D.begin(), nd.D.end());
  return nd;
}

// Empirical p-value against the pooled null [INFERRED G-08i]:
// p = (b + 1) / (N_p + 1), b = number of null scores >= observed.  [STATED: M&M ROI floor] [INFERRED G-08h]
inline std::vector<double> empirical_p(const std::vector<double>& obs, const std::vector<double>& null_sorted,
                                       bool strict = false) {
  std::vector<double> p(obs.size());
  const double Np = (double)null_sorted.size();
  for (size_t i = 0; i < obs.size(); ++i) {
    // [CODE D-12] the code counts null scores strictly greater than the observed score
    size_t b = strict ? null_sorted.end() - std::upper_bound(null_sorted.begin(), null_sorted.end(), obs[i])
                      : null_sorted.end() - std::lower_bound(null_sorted.begin(), null_sorted.end(), obs[i]);
    p[i] = ((double)b + 1.0) / (Np + 1.0);
  }
  return p;
}

}  // namespace rae
