// RAE 3.0 core -- default output set: tables, run manifest, SVG figures, HTML report.
// The same writer is used by the CLI and by the R package, so both produce identical files.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "pipeline.hpp"

namespace rae {

static const char* const RAE_VERSION = "3.0.0.9000";

inline std::string fmt(double v, const char* f = "%.6g") {
  if (v != v) return "NA";
  char b[64];
  std::snprintf(b, sizeof b, f, v);
  return b;
}
inline int64_t reg_bp_start(const Result& R, int r) { return R.G.m_start[R.ubp.start[r]]; }
inline int64_t reg_bp_end(const Result& R, int r) { return R.G.m_end[R.ubp.end[r]]; }
inline std::string arm_of(const Genome& G, int c, int64_t s, int64_t e) {
  if (G.cen_start[c] <= 0) return "";
  if (e <= G.cen_end[c]) return "p";
  if (s >= G.cen_start[c]) return "q";
  return "pq";
}

// ---- tables -------------------------------------------------------------------------------------
inline void write_regions(const Result& R, const std::string& path) {
  std::ofstream o(path);
  o << "region\tchrom\tstart\tend\tn_markers\tn_samples\tA0\tA1\tD0\tD1\terr_A0\terr_A1\terr_D0\terr_D1\t"
       "A_score\tA_err\tA_p\tA_q\tA_at_floor\tD_score\tD_err\tD_p\tD_q\tD_at_floor\t"
       "freq_gain\tfreq_amp\tfreq_loss\tfreq_homdel\tgermline_gain_normals\tgermline_loss_normals\n";
  const double floorA = 1.0 / (R.null.A.size() + 1.0), floorD = 1.0 / (R.null.D.size() + 1.0);
  for (int r = 0; r < R.ubp.size(); ++r) {
    int cnt[5] = {0, 0, 0, 0, 0};
    for (int j = 0; j < R.M.N; ++j) {
      float a0 = R.M.at(A0, r, j);
      if (a0 != a0) continue;
      cnt[discrete_state(a0, R.M.at(A1, r, j), R.M.at(D0, r, j), R.M.at(D1, r, j), R.P) + 2]++;
    }
    double n = std::max(1, R.st.n[r]);
    o << r + 1 << '\t' << R.ubp.chrom[r] + 1 << '\t' << reg_bp_start(R, r) << '\t' << reg_bp_end(R, r) << '\t'
      << R.ubp.end[r] - R.ubp.start[r] + 1 << '\t' << R.st.n[r];
    for (int k = 0; k < 4; ++k) o << '\t' << fmt(R.st.mean[k][r]);
    for (int k = 0; k < 4; ++k) o << '\t' << fmt(R.st.err[k][r]);
    o << '\t' << fmt(R.st.Ap[r]) << '\t' << fmt(R.st.errA[r]) << '\t' << fmt(R.pA[r]) << '\t' << fmt(R.qA[r])
      << '\t' << (R.pA[r] <= floorA ? 1 : 0)   // [INFERRED G-11e]
      << '\t' << fmt(R.st.Dp[r]) << '\t' << fmt(R.st.errD[r]) << '\t' << fmt(R.pD[r]) << '\t' << fmt(R.qD[r])
      << '\t' << (R.pD[r] <= floorD ? 1 : 0)
      << '\t' << fmt((cnt[3] + cnt[4]) / n) << '\t' << fmt(cnt[4] / n)      // gain incl. amp; amp  [G-12d]
      << '\t' << fmt((cnt[1] + cnt[0]) / n) << '\t' << fmt(cnt[0] / n)      // loss incl. homdel; homdel
      << '\t' << R.germ_gain[r] << '\t' << R.germ_loss[r] << '\n';
  }
}

inline void write_roi(const Result& R, const std::string& path) {
  std::ofstream o(path);
  o << "roi\tdirection\tevent\ttype\tchrom\tarm\tstart\tend\tsize_bp\tn_regions\tapex_region\tscore\terror\tq\t"
       "freq\tfreq_high\tgermline_normals_max\tcnv_n\tcnv_gain_pct\tcnv_loss_pct\tcnv_unknown_pct\tcnv_any_pct\tlikely_polymorphic";
  if (R.P.peak_independence) o << "\tresid_score\tresid_p\tindependent";   // [INFERRED G-14f]
  o << '\n';
  int id = 0;
  for (const Roi& x : R.roi) {
    int c = R.ubp.chrom[x.r_first];
    int64_t s = reg_bp_start(R, x.r_first), e = reg_bp_end(R, x.r_last);
    int hi = 0, any = 0, n = 0, germ = 0;
    for (int j = 0; j < R.M.N; ++j) {
      float a0 = R.M.at(A0, x.apex, j);
      if (a0 != a0) continue;
      ++n;
      int st = discrete_state(a0, R.M.at(A1, x.apex, j), R.M.at(D0, x.apex, j), R.M.at(D1, x.apex, j), R.P);
      if (x.direction > 0) { if (st >= 1) ++any; if (st == 2) ++hi; }
      else { if (st <= -1) ++any; if (st == -2) ++hi; }
    }
    for (int r = x.r_first; r <= x.r_last; ++r)
      germ = std::max(germ, x.direction > 0 ? R.germ_gain[r] : R.germ_loss[r]);
    o << ++id << '\t' << (x.direction > 0 ? "gain" : "loss") << '\t' << x.event + 1 << '\t'
      << (x.is_peak ? "peak" : "region") << '\t' << c + 1 << '\t' << arm_of(R.G, c, s, e) << '\t' << s << '\t' << e
      << '\t' << e - s + 1 << '\t' << x.n_regions << '\t' << x.apex + 1 << '\t' << fmt(x.score) << '\t' << fmt(x.error)
      << '\t' << fmt(x.q) << '\t' << fmt(n ? (double)any / n : NaN) << '\t' << fmt(n ? (double)hi / n : NaN)
      << '\t' << germ;
    // population CNV catalogue: reported, never used to exclude  [STATED: M&M]; flag per Taylor 2010
    CnvOverlap ov = cnv_overlap(R.G, c, s, e);
    bool poly = ov.any > R.P.cnv_polymorphic_cover && ov.n >= R.P.cnv_polymorphic_n;
    if (R.G.cnv.empty()) o << "\tNA\tNA\tNA\tNA\tNA\tNA";
    else o << '\t' << ov.n << '\t' << fmt(100 * ov.gain, "%.1f") << '\t' << fmt(100 * ov.loss, "%.1f") << '\t'
           << fmt(100 * ov.unknown, "%.1f") << '\t' << fmt(100 * ov.any, "%.1f") << '\t' << (poly ? 1 : 0);
    if (R.P.peak_independence)   // [INFERRED G-14f]
      o << '\t' << fmt(x.resid_score) << '\t' << fmt(x.resid_p) << '\t' << (x.independent < 0 ? "NA" : x.independent ? "1" : "0");
    o << '\n';
  }
}

inline std::string flag_string(uint32_t f) {
  static const char* names[] = {"sd_from_breakpoints", "dn_from_sd", "q_from_pseudo", "A0_fallback", "D0_fallback",
                                "no_halfmax", "D1_pct01", "D1_bound", "A1_bound", "ubp_excluded", "sd_cohort_median", "step_fallback", "ploidy_ambiguous"};
  std::string s;
  for (int i = 0; i < 13; ++i) if (f & (1u << i)) { if (!s.empty()) s += ","; s += names[i]; }
  return s.empty() ? "-" : s;
}

inline void write_samples(const Result& R, const std::string& path) {
  std::ofstream o(path);
  // dn = derivative noise of the paper (Eq. 1) = MAPD, median absolute difference of adjacent markers
  o << "ID\tn_segments\tin_ubp\tshift\tdn\tsd\tbandwidth\tH_gain\tH_loss\tpeak_gain\tpeak_loss\t"
       "E_A0\tbeta_A0\tE_A1\tbeta_A1\tE_D0\tbeta_D0\tE_D1\tbeta_D1\trae_fga\tfga\tn_altered\t"
       "copy_step_source\tcopy_step\tpurity\tploidy\tflags";
  if (R.P.ploidy_prior) o << "\tploidy_hard\tploidy_alt\tploidy_prior_ratio";   // [INFERRED G-14a]
  if (R.P.focal_only) o << "\tfocal_excluded_frac";                                // [INFERRED G-14d]
  o << '\n';
  for (size_t j = 0; j < R.S.size(); ++j) {
    const Sample& s = R.S[j];
    o << s.id << '\t' << s.segs.size() << '\t' << (int)R.member[j] << '\t' << fmt(s.shift) << '\t' << fmt(s.dn)
      << '\t' << fmt(s.sd_est) << '\t' << fmt(s.bw) << '\t' << fmt(s.Ha) << '\t' << fmt(s.Hd) << '\t' << fmt(s.Xa)
      << '\t' << fmt(s.Xd);
    for (int k = 0; k < 4; ++k) o << '\t' << fmt(s.E[k]) << '\t' << fmt(s.B[k]);
    static const char* src[] = {"-", "comb", "input"};
    o << '\t' << fmt(s.rae_fga) << '\t' << fmt(s.fga) << '\t' << s.n_altered << '\t' << src[s.step_src] << '\t' << fmt(s.step) << '\t' << fmt(s.purity) << '\t'
      << fmt(s.ploidy) << '\t' << flag_string(s.flags);
    if (R.P.ploidy_prior) o << '\t' << fmt(s.ploidy_hard) << '\t' << fmt(s.ploidy_alt) << '\t' << fmt(s.ploidy_prior_ratio);
    if (R.P.focal_only) o << '\t' << fmt(s.focal_excluded_frac);
    o << '\n';
  }
}

// Region x sample discrete states (Taylor 2010 rules; NOT RAE 2008): -2, -1, 0, 1, 2, NA.
inline void write_calls(const Result& R, const std::string& path) {
  std::ofstream o(path);
  o << "region\tchrom\tstart\tend";
  for (const Sample& s : R.S) o << '\t' << s.id;
  o << '\n';
  for (int r = 0; r < R.ubp.size(); ++r) {
    o << r + 1 << '\t' << R.ubp.chrom[r] + 1 << '\t' << reg_bp_start(R, r) << '\t' << reg_bp_end(R, r);
    for (int j = 0; j < R.M.N; ++j) {
      float a0 = R.M.at(A0, r, j);
      if (a0 != a0) o << "\tNA";
      else o << '\t' << discrete_state(a0, R.M.at(A1, r, j), R.M.at(D0, r, j), R.M.at(D1, r, j), R.P);
    }
    o << '\n';
  }
}

// Each sample's own segments with the four scores and the discrete state.
inline void write_segments(const Result& R, const std::string& path) {
  std::ofstream o(path);
  o << "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean.centred\tseg.sd\tA0\tA1\tD0\tD1\tstate\n";
  for (size_t j = 0; j < R.S.size(); ++j) {
    const Sample& s = R.S[j];
    for (size_t i = 0; i < s.segs.size(); ++i) {
      const Segment& g = s.segs[i];
      const SegScores& sc = R.seg_scores[j];
      o << s.id << '\t' << g.chrom + 1 << '\t' << g.bp_start << '\t' << g.bp_end << '\t' << g.nmark << '\t'
        << fmt(g.mean) << '\t' << fmt(g.sd);
      for (int k = 0; k < 4; ++k) o << '\t' << fmt(sc.s[k][i], "%.4g");
      o << '\t' << discrete_state(sc.s[A0][i], sc.s[A1][i], sc.s[D0][i], sc.s[D1][i], R.P) << '\n';
    }
  }
}

// Arm-level cohort test (arm_test only): one row per eligible arm.  [INFERRED G-14e]
inline void write_arms(const Result& R, const std::string& path) {
  std::ofstream o(path);
  o << "arm\tchrom\tarm_label\tstart\tend\tlength_bp\tn_markers\tn_samples\tA0\tA1\tD0\tD1\tscore_gain\tscore_loss\t"
       "p_gain\tq_gain\tfwer_gain\tp_loss\tq_loss\tfwer_loss\tfrac_gain\tfrac_loss\tbaseline_coverage\n";
  int id = 0;
  for (const ArmResult& r : R.arm_results) {
    const Arm& a = R.arms[r.arm];
    o << ++id << '\t' << a.chrom + 1 << '\t' << a.chrom + 1;
    if (a.arm != 'w') o << a.arm;                 // whole chromosome (centromere unknown): the number alone
    o << '\t' << a.bp_start << '\t' << a.bp_end << '\t' << a.bp_end - a.bp_start + 1 << '\t' << a.m_last - a.m_first + 1
      << '\t' << r.n_samples;
    for (int k = 0; k < 4; ++k) o << '\t' << fmt(r.mean[k]);
    o << '\t' << fmt(r.score_gain) << '\t' << fmt(r.score_loss) << '\t' << fmt(r.p_gain) << '\t' << fmt(r.q_gain) << '\t'
      << fmt(r.fwer_gain) << '\t' << fmt(r.p_loss) << '\t' << fmt(r.q_loss) << '\t' << fmt(r.fwer_loss) << '\t'
      << fmt(r.frac_gain) << '\t' << fmt(r.frac_loss) << '\t' << fmt(r.baseline_coverage) << '\n';
  }
}

inline std::string json_escape(const std::string& s) {
  std::string o;
  for (char c : s) { if (c == '"' || c == '\\') o += '\\'; if (c == '\n') { o += "\\n"; continue; } o += c; }
  return o;
}

// Run manifest: every parameter, the seed, and which inferred fallbacks fired.
inline void write_manifest(const Result& R, const std::string& path, const std::string& input) {
  std::ofstream o(path);
  const Params& P = R.P;
  int nf[13] = {0};
  for (const Sample& s : R.S) for (int i = 0; i < 13; ++i) if (s.flags & (1u << i)) nf[i]++;
  int nsigA = 0, nsigD = 0, nreg = 0, npeak = 0;
  for (int r = 0; r < R.ubp.size(); ++r) { nsigA += R.qA[r] <= P.fdr; nsigD += R.qD[r] <= P.fdr; }
  for (const Roi& x : R.roi) { if (x.is_peak) ++npeak; else ++nreg; }
  o << "{\n  \"tool\": \"rae\",\n  \"version\": \"" << RAE_VERSION << "\",\n"
    << "  \"method\": \"Reimplementation of RAE, Taylor et al. 2008 PLoS ONE 3(9):e3179; see docs/spec/02_gap_register.md for every inferred choice\",\n"
    << "  \"input\": \"" << json_escape(input) << "\",\n"
    << "  \"genome\": {\"build\": \"" << R.G.build << "\", \"n_markers\": " << R.G.n_markers()
    << ", \"virtual_grid\": " << (R.G.virtual_grid ? "true" : "false") << ", \"grid_spacing\": " << json_num(R.G.grid_spacing)
    << ", \"hotspot_table\": " << (R.G.has_hotspots ? "true" : "false")
    << ", \"cnv_catalogue\": \"" << json_escape(R.G.cnv_source) << "\", \"cnv_catalogue_rows\": " << R.G.cnv.size() << "},\n"
    << "  \"input_log\": {\"rows_read\": " << R.log.rows_read << ", \"rows_non_autosomal_dropped\": " << R.log.rows_non_autosomal
    << ", \"rows_unparseable\": " << R.log.rows_bad << ", \"segments_unmappable\": " << R.log.segs_dropped_unmappable
    << ", \"segments_duplicate_start\": " << R.log.segs_dropped_duplicate << ", \"seg_sd_column\": " << (R.log.has_sd ? "true" : "false") << "},\n"
    << "  \"parameters\": {" << params_json(P) << "},\n"
    << "  \"inferred_fallbacks_fired\": {\n"
    << "    \"G-01a_sd_from_breakpoints_samples\": " << nf[0] << ", \"G-01a_sd_cohort_median_samples\": " << nf[10]
    << ", \"G-01b_dn_from_sd_samples\": " << nf[1] << ", \"G-01c_percentiles_from_pseudo_probes_samples\": " << nf[2]
    << ",\n    \"G-01d_virtual_marker_grid\": " << (R.G.virtual_grid ? "true" : "false")
    << ", \"G-08k_uniform_cleavage\": " << (R.G.has_hotspots ? "false" : "true")
    << ", \"G-03f_no_halfmax_samples\": " << nf[5] << ", \"G-13b_copy_step_fallback_samples\": " << nf[11];
  if (P.ploidy_prior) o << ", \"G-14b_ploidy_ambiguous_samples\": " << nf[12];   // [INFERRED G-14b]
  o << "\n  },\n"
    << "  \"summary\": {\"n_samples\": " << R.S.size() << ", \"n_normals\": " << R.n_normals << ", \"ubp_excluded_samples\": " << nf[9]
    << ", \"ubp_regions\": " << R.ubp.size() << ", \"median_segment_bp\": " << json_num(R.median_segment_bp)
    << ",\n    \"A0_fallback_samples\": " << nf[3] << ", \"D0_fallback_samples\": " << nf[4] << ", \"D1_pct01_samples\": " << nf[6]
    << ", \"D1_bound_samples\": " << nf[7] << ", \"A1_bound_samples\": " << nf[8]
    << ",\n    \"null_rounds\": " << R.null_info.rounds << ", \"null_scores\": " << R.null.A.size()
    << ", \"cleaved_regions_first_round\": " << (R.null_info.cleaved_regions.empty() ? 0 : R.null_info.cleaved_regions[0])
    << ", \"p_floor\": " << json_num(1.0 / (R.null.A.size() + 1.0))
    << ",\n    \"significant_regions_gain\": " << nsigA << ", \"significant_regions_loss\": " << nsigD
    << ", \"roi_events\": " << nreg << ", \"roi_peaks\": " << npeak;
  if (P.arm_test) {   // [INFERRED G-14e]
    int ag = 0, al = 0;
    for (const ArmResult& r : R.arm_results) { ag += r.q_gain <= P.fdr; al += r.q_loss <= P.fdr; }
    o << ", \"arms_significant_gain\": " << ag << ", \"arms_significant_loss\": " << al;
  }
  o << "},\n"
    << "  \"warnings\": [";
  for (size_t i = 0; i < R.log.warnings.size(); ++i) o << (i ? ", " : "") << '"' << json_escape(R.log.warnings[i]) << '"';
  o << "]\n}\n";
}

// Null-model internals for the diagnostic figures (Taylor 2008 Fig. 3B, 3C), in long format:
// item chrom value. chrom_len, cen_start, cen_end (bp; 0 = unknown); hotspot (centre, bp); sampled_cut
// (first null round, bp); ubp_bp and cleaved_bp (region sizes of the unified breakpoint profile and of the
// first hotspot-cleaved null round, chromosome of the region's first marker).
inline void write_null(const Result& R, const std::string& path) {
  std::ofstream o(path);
  const Genome& G = R.G;
  o << "item\tchrom\tvalue\n";
  for (int c = 0; c < G.n_chrom; ++c) o << "chrom_len\t" << c + 1 << '\t' << G.chrom_len[c] << '\n';
  for (int c = 0; c < G.n_chrom; ++c) o << "cen_start\t" << c + 1 << '\t' << G.cen_start[c] << '\n';
  for (int c = 0; c < G.n_chrom; ++c) o << "cen_end\t" << c + 1 << '\t' << G.cen_end[c] << '\n';
  for (int c = 0; c < G.n_chrom && c < (int)G.hotspots.size(); ++c)
    for (int64_t x : G.hotspots[c]) o << "hotspot\t" << c + 1 << '\t' << x << '\n';
  for (int32_t m : R.null_info.first_cuts)
    if (m >= 0 && m < G.n_markers()) o << "sampled_cut\t" << G.m_chrom[m] + 1 << '\t' << G.m_start[m] << '\n';
  auto sizes = [&](const char* item, const Regions& r) {
    for (int i = 0; i < r.size(); ++i)
      o << item << '\t' << G.m_chrom[r.start[i]] + 1 << '\t' << G.m_end[r.end[i]] - G.m_start[r.start[i]] + 1 << '\n';
  };
  sizes("ubp_bp", R.ubp);
  sizes("cleaved_bp", R.cleaved_example);
}

}  // namespace rae
#include "figures.hpp"

namespace rae {
// The fixed default output set. Returns the list of files written.
inline std::vector<std::string> write_outputs(const Result& R, const std::string& prefix, const std::string& input) {
  std::vector<std::string> f;
  std::filesystem::path parent = std::filesystem::path(prefix).parent_path();
  if (!parent.empty()) std::filesystem::create_directories(parent);
  auto p = [&](const char* suffix) { f.push_back(prefix + suffix); return f.back(); };
  write_regions(R, p(".regions.tsv"));
  write_roi(R, p(".roi.tsv"));
  write_samples(R, p(".samples.tsv"));
  write_segments(R, p(".segments.tsv"));
  write_calls(R, p(".calls.tsv"));
  write_null(R, p(".null.tsv"));
  if (R.P.arm_test) write_arms(R, p(".arms.tsv"));   // a ninth file, only with the arm test  [G-14e]
  else std::filesystem::remove(prefix + ".arms.tsv");  // never leave one from an earlier run under this prefix
  write_manifest(R, p(".run.json"), input);
  write_report(R, p(".report.html"), input);
  for (const std::string& x : f)  // never report a file that was not written
    if (!std::filesystem::exists(x) || std::filesystem::file_size(x) == 0) throw std::runtime_error("could not write " + x);
  return f;
}
}  // namespace rae
