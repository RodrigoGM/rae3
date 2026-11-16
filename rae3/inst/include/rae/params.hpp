// RAE 3.0 core -- one table of every parameter, used to set them by name (CLI --set, R list) and to
// write them all into the run manifest.
#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include "types.hpp"

namespace rae {

// A JSON number: shortest of %.15g..%.17g that reads back exactly; null for NaN and infinities.
inline std::string json_num(double v) {
  if (!std::isfinite(v)) return "null";
  char b[40];
  for (int p = 15; p <= 17; ++p) {
    std::snprintf(b, sizeof b, "%.*g", p, v);
    if (std::strtod(b, nullptr) == v) break;
  }
  return b;
}

#define RAE_PARAM_TABLE(X) \
  X(ubp_min_markers) X(density_n) X(density_from) X(density_to) X(density_adjust) X(peak_window_lo) \
  X(peak_window_hi) X(sigmoid_at_half) X(beta_fallback_gain) X(beta_fallback_loss) X(d1_pct_first) \
  X(d1_pct_second) X(d1_bound) X(a1_bound) X(beta_a1) X(n_perm) X(recleave_every) X(fdr) X(peak_q) \
  X(shoulder_factor) X(germline_min_score) X(germline_min_normals) X(null_positions) X(ubp_outlier_k) \
  X(split_at_centromere) X(error_sem) X(sd_from_breakpoints_c) X(marker_sd) X(grid_spacing) \
  X(uniform_cleave_spacing) X(density_bw) X(density_noise_sem) X(d0_midpoint_zero) X(beta_min) \
  X(halfmax_min) X(peak_span) X(diploid_nearest_zero) X(centre_window) X(e_inflate) X(ubp_outlier_log) \
  X(hotspot_kde) X(cleave_code) X(p_strict) X(roi_code) X(roi_error_factor) X(region_weight_code) \
  X(call_threshold) X(amp_threshold) X(burden_threshold) X(fga_threshold) X(fga_threshold_loss) X(fga_assembly) X(cnv_polymorphic_cover) X(cnv_polymorphic_n) \
  X(copy_step) X(comb_step_min) X(comb_step_max) X(comb_step_by) X(comb_rel) X(comb_tooth_min_frac) \
  X(amp_cn_low) X(amp_cn_high) X(amp_ploidy_split) X(homdel_cn) X(a1_logistic) \
  X(ploidy_prior) X(ploidy_prior_wgd) X(ploidy_prior_mu1) X(ploidy_prior_mu2) X(ploidy_prior_sd) X(ploidy_prior_ambig) \
  X(arm_baseline) X(arm_min_markers) X(focal_only) X(focal_max_bp) X(arm_test) X(arm_n_perm) \
  X(peak_independence) X(peak_run_min)

inline void assign_param(int& f, double v) { f = (int)v; }
inline void assign_param(bool& f, double v) { f = v != 0; }
inline void assign_param(double& f, double v) { f = v; }

// Returns false when the name is unknown.
inline bool set_param(Params& P, const std::string& name, double value) {
  if (name == "seed") { P.seed = (uint64_t)value; return true; }
#define X(f) if (name == #f) { assign_param(P.f, value); return true; }
  RAE_PARAM_TABLE(X)
#undef X
  return false;
}

inline bool set_preset(Params& P, const std::string& name) {
  // A preset resets every parameter, so it must come before individual overrides.
  if (name == "sequencing") { P = Params(); return true; }
  if (name == "taylor2008") { P = Params(); P.apply_taylor2008(); return true; }
  if (name == "cgh") { P = Params(); P.apply_cgh(); return true; }
  return false;
}

inline std::string params_json(const Params& P) {
  std::ostringstream o;
  o << "\"preset\": \"" << P.preset << "\", \"seed\": " << P.seed;
  auto num = [&](const char* n, double v) { o << ", \"" << n << "\": " << json_num(v); };
#define X(f) num(#f, (double)P.f);
  RAE_PARAM_TABLE(X)
#undef X
  return o.str();
}

}  // namespace rae
