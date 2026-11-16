// Rcpp bindings to the header-only RAE core (inst/include/rae). The whole-pipeline entry point
// writes the same default output set as the CLI; the step functions expose each module so it can be
// validated on its own.
#include <Rcpp.h>
#include "rae/arms.hpp"
#include "rae/report.hpp"

using namespace Rcpp;

static rae::Params params_from_list(const List& L) {
  rae::Params P;
  if (L.containsElementNamed("preset") && !Rf_isNull(L["preset"])) {
    std::string pr = as<std::string>(L["preset"]);
    if (!rae::set_preset(P, pr)) stop("preset must be \"sequencing\", \"cgh\" or \"taylor2008\"");
  }
  CharacterVector nm = L.names();
  for (R_xlen_t i = 0; i < L.size(); ++i) {
    std::string name = as<std::string>(nm[i]);
    if (name == "preset" || Rf_isNull(L[i])) continue;
    if (!rae::set_param(P, name, as<double>(L[i]))) stop("unknown RAE parameter: %s", name);
  }
  return P;
}

static std::vector<rae::SegRow> rows_from_df(const DataFrame& d, rae::InputLog& log) {
  CharacterVector id = d[0], ch = d[1];
  NumericVector st = d[2], en = d[3], nm = d[4], mean = d[5];
  bool has_sd = d.size() >= 7;
  NumericVector sd = has_sd ? as<NumericVector>(d[6]) : NumericVector(0);
  std::vector<rae::SegRow> rows;
  for (R_xlen_t i = 0; i < id.size(); ++i) {
    log.rows_read++;
    int c = rae::parse_chrom(as<std::string>(ch[i]));
    if (c < 0) { log.rows_non_autosomal++; continue; }
    if (NumericVector::is_na(st[i]) || NumericVector::is_na(en[i]) || NumericVector::is_na(mean[i]) || !std::isfinite(mean[i]) ||
        !rae::row_in_range(st[i], en[i], NumericVector::is_na(nm[i]) ? rae::NaN : (double)nm[i])) { log.rows_bad++; continue; }
    double s = has_sd && !NumericVector::is_na(sd[i]) && sd[i] >= 0 && std::isfinite(sd[i]) ? sd[i] : rae::NaN;
    if (s == s) log.has_sd = true;
    rows.push_back({as<std::string>(id[i]), c, (int64_t)st[i], (int64_t)en[i],
                    NumericVector::is_na(nm[i]) ? 1 : (int)nm[i], mean[i], s});
  }
  return rows;
}

// bins = chrom (character), start, end; NULL gives no bins (a virtual grid is then reconstructed).
static std::vector<rae::BinRow> bins_from_df(const Nullable<DataFrame>& bins) {
  std::vector<rae::BinRow> b;
  if (bins.isNull()) return b;
  DataFrame bd(bins);
  CharacterVector ch = bd[0];
  NumericVector st = bd[1], en = bd[2];
  for (R_xlen_t i = 0; i < ch.size(); ++i) {
    int c = rae::parse_chrom(as<std::string>(ch[i]));
    if (c >= 0 && !NumericVector::is_na(st[i]) && !NumericVector::is_na(en[i]) && rae::row_in_range(st[i], en[i], rae::NaN))
      b.push_back({c, (int64_t)st[i], (int64_t)en[i]});
  }
  return b;
}

static std::map<std::string, rae::StatRow> stats_from_df(const Nullable<DataFrame>& n) {
  std::map<std::string, rae::StatRow> m;
  if (n.isNull()) return m;
  DataFrame d(n);
  CharacterVector id = d["ID"];
  auto col = [&](const char* n) { return d.containsElementNamed(n) ? as<NumericVector>(d[n]) : NumericVector(id.size(), NA_REAL); };
  NumericVector dn = col("dn"), q1 = col("q01"), q25 = col("q025"), pu = col("purity"), pl = col("ploidy");
  for (R_xlen_t i = 0; i < id.size(); ++i) {
    rae::StatRow s;
    s.dn = NumericVector::is_na(dn[i]) ? rae::NaN : dn[i];
    s.q01 = NumericVector::is_na(q1[i]) ? rae::NaN : q1[i];
    s.q025 = NumericVector::is_na(q25[i]) ? rae::NaN : q25[i];
    s.purity = NumericVector::is_na(pu[i]) ? rae::NaN : pu[i];
    s.ploidy = NumericVector::is_na(pl[i]) ? rae::NaN : pl[i];
    m[as<std::string>(id[i])] = s;
  }
  return m;
}

// Whole pipeline. Data frames are positional: seg = ID, chrom (character), start, end, num.mark,
// seg.mean [, seg.sd]; bins = chrom (character), start, end. Writes the default output set to `prefix`.
// Returns the file names and, for the R plotting functions only, the null-model example of Taylor 2008
// Fig. 3B/C: hotspot centres, the cut positions of the first randomisation, and region sizes (bp) of the
// unified breakpoint profile and of its first hotspot-cleaved counterpart. The command-line output is unchanged.
// [[Rcpp::export]]
List rae_run_cpp(DataFrame seg, Nullable<DataFrame> bins, Nullable<DataFrame> stats,
                            Nullable<DataFrame> normals, Nullable<DataFrame> normal_stats,
                            std::string build, List params, std::string prefix, std::string input_label,
                            std::string cnv = "default") {
  rae::Params P = params_from_list(params);
  rae::InputLog log, nlog;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  std::vector<rae::BinRow> b = bins_from_df(bins);
  std::vector<rae::SegRow> nrows;
  if (normals.isNotNull()) nrows = rows_from_df(DataFrame(normals), nlog);
  std::map<std::string, rae::StatRow> st = stats_from_df(stats);
  rae::Result R = rae::run(rows, b, st, build, P, log, nrows, stats_from_df(normal_stats), cnv);
  std::vector<std::string> files = rae::write_outputs(R, prefix, input_label);
  // the null-model internals are in <prefix>.null.tsv, read back by rae_read()
  return List::create(_["files"] = wrap(files));
}

// Test hook (not exported): runs the pipeline and summarises the pooled null scores (sorted ascending): max and 99th
// percentile for gain (A) and loss (D).
// [[Rcpp::export]]
NumericVector rae_null_summary_cpp(DataFrame seg, std::string build, List params) {
  rae::Params P = params_from_list(params);
  rae::InputLog log;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  rae::Result R = rae::run(rows, {}, {}, build, P, log);
  auto q = [](const std::vector<double>& v, double p) { return v.empty() ? NA_REAL : rae::quantile7_sorted(v, p); };
  return NumericVector::create(_["A_max"] = R.null.A.empty() ? NA_REAL : R.null.A.back(), _["A_q99"] = q(R.null.A, 0.99),
                               _["D_max"] = R.null.D.empty() ? NA_REAL : R.null.D.back(), _["D_q99"] = q(R.null.D, 0.99));
}

// ---- step functions -------------------------------------------------------------------------------

// [[Rcpp::export]]
double rae_bw_nrd0(NumericVector x) { return rae::bw_nrd0(as<std::vector<double>>(x)); }

// [[Rcpp::export]]
List rae_density(NumericVector x, double bw, double from = -1, double to = 1, int n = 2048) {
  rae::Density d = rae::density_r(as<std::vector<double>>(x), bw, from, to, n);
  return List::create(_["x"] = d.x, _["y"] = d.y, _["bw"] = d.bw);
}

// component: "A0", "A1", "D0", "D1". Eq. 2 of the paper.
// [[Rcpp::export]]
NumericVector rae_sigmoid(NumericVector a, double E, double beta, std::string component = "A0") {
  NumericVector out(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i)
    out[i] = component == "A1" ? rae::sigmoid_a1(a[i], E, beta) : rae::sigmoid(a[i], E, beta);
  return out;
}

// First derivative used in the analytical error (Methods S1-C).
// [[Rcpp::export]]
NumericVector rae_dsigmoid(NumericVector a, double E, double beta, std::string component = "A0") {
  NumericVector out(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i)
    out[i] = component == "A1" ? (a[i] > E ? 2 * rae::dsigmoid(a[i], E, beta) : 0.0) : rae::dsigmoid(a[i], E, beta);
  return out;
}

// Parameterise one sample from its segments. dn, q01, q025 may be NA (then the documented fallbacks fire).
// [[Rcpp::export]]
List rae_parameterize_cpp(NumericVector mean, NumericVector sd, IntegerVector nmark, double dn, double q01,
                          double q025, List params, double stream = 0) {
  rae::Params P = params_from_list(params);
  rae::Sample s;
  s.id = "sample";
  for (R_xlen_t i = 0; i < mean.size(); ++i) {
    rae::Segment g;
    g.start = (int32_t)i; g.end = (int32_t)i; g.chrom = 0; g.nmark = std::max(1, nmark[i]);
    g.mean = mean[i]; g.sd = NumericVector::is_na(sd[i]) ? rae::NaN : sd[i];
    s.segs.push_back(g);
  }
  s.dn = NumericVector::is_na(dn) ? rae::NaN : dn;
  s.q01 = NumericVector::is_na(q01) ? rae::NaN : q01;
  s.q025 = NumericVector::is_na(q025) ? rae::NaN : q025;
  rae::resolve_sd(s, P, 0.1);
  rae::parameterize_sample(s, P, (uint64_t)stream);
  NumericVector E(4), B(4);
  for (int k = 0; k < 4; ++k) { E[k] = s.E[k]; B[k] = s.B[k]; }
  E.names() = CharacterVector::create("A0", "A1", "D0", "D1");
  B.names() = CharacterVector::create("A0", "A1", "D0", "D1");
  return List::create(_["E"] = E, _["beta"] = B, _["shift"] = s.shift, _["dn"] = s.dn, _["sd"] = s.sd_est,
                      _["bw"] = s.bw, _["H_gain"] = s.Ha, _["H_loss"] = s.Hd, _["peak_gain"] = s.Xa,
                      _["peak_loss"] = s.Xd, _["flags"] = rae::flag_string(s.flags),
                      _["density_x"] = s.dens_x, _["density_y"] = s.dens_y);
}

// Unified breakpoint profile from per-sample segment starts on a single chromosome of n_markers
// markers (1-based marker indices). `starts` and `nmark` are lists with one vector per sample.
// [[Rcpp::export]]
DataFrame rae_ubp_cpp(List starts, List nmark, int n_markers, List params) {
  rae::Params P = params_from_list(params);
  rae::Genome G;
  G.n_chrom = 1;
  G.chrom_len.assign(1, n_markers); G.cen_start.assign(1, 0); G.cen_end.assign(1, 0);
  G.chrom_first = {0, n_markers};
  for (int m = 0; m < n_markers; ++m) { G.m_chrom.push_back(0); G.m_start.push_back(m + 1); G.m_end.push_back(m + 1); }
  std::vector<rae::Profile> prof(starts.size());
  for (R_xlen_t j = 0; j < starts.size(); ++j) {
    IntegerVector s = starts[j], n = nmark[j];
    for (R_xlen_t i = 0; i < s.size(); ++i) {
      int end = i + 1 < s.size() ? s[i + 1] - 2 : n_markers - 1;
      prof[j].push_back({s[i] - 1, end, (int32_t)i, n[i] >= P.ubp_min_markers});
    }
  }
  std::vector<char> member(starts.size(), 1);
  rae::Regions reg = rae::build_ubp(prof, member, G, P);
  IntegerVector a(reg.size()), b(reg.size());
  for (int r = 0; r < reg.size(); ++r) { a[r] = reg.start[r] + 1; b[r] = reg.end[r] + 1; }
  return DataFrame::create(_["start"] = a, _["end"] = b);
}

// [[Rcpp::export]]
NumericVector rae_bh(NumericVector p) { return wrap(rae::bh_adjust(as<std::vector<double>>(p))); }

// [[Rcpp::export]]
NumericVector rae_empirical_p(NumericVector observed, NumericVector null_scores) {
  std::vector<double> nl = as<std::vector<double>>(null_scores);
  std::sort(nl.begin(), nl.end());
  return wrap(rae::empirical_p(as<std::vector<double>>(observed), nl));
}

// Regions of interest from per-region score, error, q and chromosome (one direction).
// [[Rcpp::export]]
DataFrame rae_roi_cpp(NumericVector score, NumericVector error, NumericVector q, IntegerVector chrom, List params) {
  rae::Params P = params_from_list(params);
  std::vector<rae::Roi> roi = rae::find_roi(as<std::vector<double>>(score), as<std::vector<double>>(error),
                                            as<std::vector<double>>(q), as<std::vector<int32_t>>(chrom), +1, P);
  size_t n = roi.size();
  IntegerVector ev(n), first(n), last(n), apex(n);
  LogicalVector peak(n);
  for (size_t i = 0; i < n; ++i) {
    ev[i] = roi[i].event + 1; first[i] = roi[i].r_first + 1; last[i] = roi[i].r_last + 1;
    apex[i] = roi[i].apex + 1; peak[i] = roi[i].is_peak;
  }
  return DataFrame::create(_["event"] = ev, _["is_peak"] = peak, _["first"] = first, _["last"] = last, _["apex"] = apex);
}

// Uniform draws from the core generator (for reproducibility tests).
// [[Rcpp::export]]
NumericVector rae_rng_unif(int n, double seed, double stream = 0) {
  rae::Rng rng((uint64_t)seed, (uint64_t)stream);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = rng.unif();
  return out;
}

// Bundled genome builds.
// [[Rcpp::export]]
CharacterVector rae_builds() {
  CharacterVector out(rae::annot::n_builds);
  for (int i = 0; i < rae::annot::n_builds; ++i) out[i] = rae::annot::builds[i].name;
  return out;
}

// ---- arm geometry and baselines (G-14c). Internal test hooks; not exported from NAMESPACE. ---------------

static rae::Genome arm_genome(const DataFrame& bins, const std::string& build, const rae::Params& P,
                              const std::vector<rae::SegRow>& rows) {
  std::vector<rae::BinRow> b = bins_from_df(Nullable<DataFrame>(bins));
  if (b.empty()) stop("bins must hold at least one marker");
  rae::InputLog log;
  return rae::make_genome(build, rows, b, P, log);
}

// Arms of a bundled build on a marker grid `bins` (chrom, start, end). Marker indices are 1-based and global.
// [[Rcpp::export]]
DataFrame rae_arm_layout_cpp(std::string build, int min_markers, DataFrame bins) {
  rae::Params P;
  P.arm_min_markers = min_markers;
  rae::Genome G = arm_genome(bins, build, P, {});
  std::vector<rae::Arm> arms = rae::arm_layout(G, P);
  size_t n = arms.size();
  IntegerVector chrom(n), mf(n), ml(n), nm(n);
  NumericVector bs(n), be(n);
  LogicalVector el(n);
  CharacterVector arm(n);
  for (size_t i = 0; i < n; ++i) {
    chrom[i] = arms[i].chrom + 1; arm[i] = std::string(1, arms[i].arm);
    mf[i] = arms[i].m_first + 1; ml[i] = arms[i].m_last + 1; nm[i] = arms[i].m_last - arms[i].m_first + 1;
    bs[i] = (double)arms[i].bp_start; be[i] = (double)arms[i].bp_end; el[i] = arms[i].eligible;
  }
  return DataFrame::create(_["chrom"] = chrom, _["arm"] = arm, _["m_first"] = mf, _["m_last"] = ml, _["n_markers"] = nm,
                           _["bp_start"] = bs, _["bp_end"] = be, _["eligible"] = el, _["stringsAsFactors"] = false);
}

// Segments of one cohort cut at the centromere; returns the seg table (ID, chrom, loc.start, loc.end, num.mark, seg.mean, seg.sd).
// [[Rcpp::export]]
DataFrame rae_split_segments_cpp(DataFrame seg, DataFrame bins, std::string build, int min_markers) {
  rae::Params P;
  P.arm_min_markers = min_markers;
  rae::InputLog log;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  rae::Genome G = arm_genome(bins, build, P, rows);
  std::vector<rae::Sample> S = rae::make_samples(rows, G, {}, log);
  rae::split_segments_at_centromeres(S, G, rae::arm_layout(G, P));
  std::vector<std::string> id;
  std::vector<int> chrom, nm;
  std::vector<double> st, en, mean, sd;
  for (const rae::Sample& s : S)
    for (const rae::Segment& g : s.segs) {
      id.push_back(s.id); chrom.push_back(g.chrom + 1); st.push_back((double)g.bp_start); en.push_back((double)g.bp_end);
      nm.push_back(g.nmark); mean.push_back(g.mean); sd.push_back(g.sd);
    }
  return DataFrame::create(_["ID"] = id, _["chrom"] = chrom, _["loc.start"] = st, _["loc.end"] = en, _["num.mark"] = nm,
                           _["seg.mean"] = mean, _["seg.sd"] = sd, _["stringsAsFactors"] = false);
}

// Per-arm level and coverage of one sample's segments: cut at the centromere, then arm_baselines(). A ladder is
// given by step, lad_m, lad_r0 (copy units; step NA: log2 scale); E and B (length 4: A0, A1, D0, D1) are the detectors
// used for coverage on the log2 scale.
// [[Rcpp::export]]
DataFrame rae_arm_baselines_cpp(DataFrame seg, DataFrame bins, std::string build, List params, double step,
                                double lad_m, double lad_r0, NumericVector E, NumericVector B) {
  rae::Params P = params_from_list(params);
  rae::InputLog log;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  rae::Genome G = arm_genome(bins, build, P, rows);
  std::vector<rae::Sample> S = rae::make_samples(rows, G, {}, log);
  if (S.empty()) stop("no segments");
  std::vector<rae::Arm> arms = rae::arm_layout(G, P);
  rae::split_segments_at_centromeres(S, G, arms);
  rae::Sample& s = S[0];
  if (!std::isnan(step)) { s.step_src = 2; s.step = step; s.lad_m = lad_m; s.lad_r0 = lad_r0; }
  if (E.size() == 4 && B.size() == 4) for (int k = 0; k < 4; ++k) { s.E[k] = E[k]; s.B[k] = B[k]; }
  s.a1_logistic = true;
  std::vector<rae::ArmLevel> lv = rae::arm_baselines(s, arms, P);
  std::vector<int> chrom;
  std::vector<std::string> arm;
  std::vector<double> level, cov;
  for (size_t i = 0; i < arms.size(); ++i) {
    if (std::isnan(lv[i].level)) continue;
    chrom.push_back(arms[i].chrom + 1); arm.push_back(std::string(1, arms[i].arm));
    level.push_back(lv[i].level); cov.push_back(lv[i].coverage);
  }
  return DataFrame::create(_["chrom"] = chrom, _["arm"] = arm, _["level"] = level, _["coverage"] = cov,
                           _["stringsAsFactors"] = false);
}

// Focal scores of one sample (G-14d): cut at the centromere, arm baselines, then focal_scores(). The ladder is given by
// step, lad_m, lad_r0 and ploidy (step NA: log2 scale); fallback sets F_STEP_FALLBACK, which forces the log2 scale;
// E and B (A0, A1, D0, D1) are the detectors. Returns the split segments with their scores; attribute
// "excluded_frac" holds the excluded share of the covered span.
// [[Rcpp::export]]
DataFrame rae_focal_scores_cpp(DataFrame seg, DataFrame bins, std::string build, List params, double step, double lad_m,
                               double lad_r0, double ploidy, bool fallback, NumericVector E, NumericVector B) {
  rae::Params P = params_from_list(params);
  rae::InputLog log;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  rae::Genome G = arm_genome(bins, build, P, rows);
  std::vector<rae::Sample> S = rae::make_samples(rows, G, {}, log);
  if (S.empty()) stop("no segments");
  if (E.size() != 4 || B.size() != 4) stop("E and B need four values (A0, A1, D0, D1)");
  std::vector<rae::Arm> arms = rae::arm_layout(G, P);
  rae::split_segments_at_centromeres(S, G, arms);
  rae::Sample& s = S[0];
  if (!std::isnan(step)) { s.step_src = 2; s.step = step; s.lad_m = lad_m; s.lad_r0 = lad_r0; s.ploidy = ploidy; }
  if (fallback) s.flags |= rae::F_STEP_FALLBACK;
  for (int k = 0; k < 4; ++k) { s.E[k] = E[k]; s.B[k] = B[k]; }
  s.a1_logistic = true;
  s.sd_est = 0.1;
  double ex = rae::NaN;
  rae::SegScores sc = rae::focal_scores(s, arms, rae::arm_baselines(s, arms, P), P, &ex);
  std::vector<int> chrom, nm;
  std::vector<double> st, en, mean;
  for (const rae::Segment& g : s.segs) {
    chrom.push_back(g.chrom + 1); st.push_back((double)g.bp_start); en.push_back((double)g.bp_end);
    nm.push_back(g.nmark); mean.push_back(g.mean);
  }
  auto col = [&](int k) { return std::vector<double>(sc.s[k].begin(), sc.s[k].end()); };
  DataFrame out = DataFrame::create(_["chrom"] = chrom, _["start"] = st, _["end"] = en, _["nmark"] = nm, _["mean"] = mean,
                                    _["A0"] = col(rae::A0), _["A1"] = col(rae::A1), _["D0"] = col(rae::D0),
                                    _["D1"] = col(rae::D1), _["stringsAsFactors"] = false);
  out.attr("excluded_frac") = ex;
  return out;
}

// The estimators alone, on per-segment values, weights (markers), chromosome and arm label ("p", "q", "w"), in
// positional order; method 0 weighted median, 1 two-level fit. One row per (chrom, arm) present, in order of
// first appearance. Coverage is on the copy-unit scale: the weight share within 0.5 of the arm's level.
// [[Rcpp::export]]
DataFrame rae_arm_baseline_cpp(NumericVector values, NumericVector weights, IntegerVector chrom_id, CharacterVector arm_label,
                               int method) {
  const R_xlen_t n = values.size();
  if (weights.size() != n || chrom_id.size() != n || arm_label.size() != n) stop("values, weights, chrom and arm must have the same length");
  std::map<std::pair<int, std::string>, int> idx;
  std::vector<std::pair<int, std::string>> keys;
  std::vector<int> chrom(n), arm(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    std::pair<int, std::string> k(chrom_id[i], as<std::string>(arm_label[i]));
    auto it = idx.find(k);
    if (it == idx.end()) { it = idx.emplace(k, (int)keys.size()).first; keys.push_back(k); }
    chrom[i] = chrom_id[i]; arm[i] = it->second;
  }
  std::vector<double> v = as<std::vector<double>>(values), w = as<std::vector<double>>(weights);
  std::vector<double> lev = rae::arm_levels(chrom, arm, v, w, (int)keys.size(), method);
  std::vector<double> at(keys.size(), 0.0), tot(keys.size(), 0.0);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!std::isfinite(v[i]) || !(w[i] > 0)) continue;
    tot[arm[i]] += w[i];
    if (std::fabs(v[i] - lev[arm[i]]) < 0.5) at[arm[i]] += w[i];
  }
  IntegerVector ch(keys.size());
  CharacterVector al(keys.size());
  NumericVector lv(keys.size()), cov(keys.size());
  for (size_t a = 0; a < keys.size(); ++a) {
    ch[a] = keys[a].first; al[a] = keys[a].second; lv[a] = lev[a]; cov[a] = tot[a] > 0 ? at[a] / tot[a] : NA_REAL;
  }
  return DataFrame::create(_["chrom"] = ch, _["arm"] = al, _["level"] = lv, _["coverage"] = cov, _["stringsAsFactors"] = false);
}

// Per-sample, per-arm state of a whole run (G-14c, G-14e), for validation of the arm baseline: runs the pipeline
// (rae::run) with `params`, which must switch on arm_test or focal_only, and returns one row per sample and arm of the
// layout with the arm baseline, its coverage and the four arm scores the arm test averages (arm_scores(); NaN where
// the sample has no baseline on the arm). stats = ID, dn, q01, q025 [, purity, ploidy], or NULL. Writes nothing.
// [[Rcpp::export]]
DataFrame rae_arm_states_cpp(DataFrame seg, DataFrame bins, Nullable<DataFrame> stats, std::string build, List params) {
  rae::Params P = params_from_list(params);
  if (!P.arm_test && !P.focal_only) stop("rae_arm_states_cpp needs arm_test or focal_only");
  rae::InputLog log;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  rae::Result R = rae::run(rows, bins_from_df(Nullable<DataFrame>(bins)), stats_from_df(stats), build, P, log);
  std::vector<std::string> id, arm;
  std::vector<int> chrom;
  std::vector<bool> el;
  std::vector<double> level, cov, sc[4];
  for (size_t j = 0; j < R.S.size(); ++j)
    for (size_t a = 0; a < R.arms.size(); ++a) {
      const bool has = a < R.baselines[j].size();
      const double b = has ? R.baselines[j][a].level : rae::NaN;
      double x[4];
      rae::arm_scores(R.S[j], b, x);
      id.push_back(R.S[j].id); chrom.push_back(R.arms[a].chrom + 1); arm.push_back(std::string(1, R.arms[a].arm));
      el.push_back(R.arms[a].eligible); level.push_back(b); cov.push_back(has ? R.baselines[j][a].coverage : rae::NaN);
      for (int k = 0; k < 4; ++k) sc[k].push_back(x[k]);
    }
  return DataFrame::create(_["ID"] = id, _["chrom"] = chrom, _["arm"] = arm, _["eligible"] = el, _["level"] = level,
                           _["coverage"] = cov, _["A0"] = sc[rae::A0], _["A1"] = sc[rae::A1], _["D0"] = sc[rae::D0],
                           _["D1"] = sc[rae::D1], _["stringsAsFactors"] = false);
}

// Test hook for the masked recompute of the peak-independence test (G-14f): runs the pipeline (rae::run) with
// `params` and recomputes every sample's region means through rae::masked_region_means with nothing masked. Returns
// list(observed, recomputed): two (samples x regions) x 4 matrices (A0, A1, D0, D1) of R.M and the recompute, which
// must be identical. bins and stats (ID, dn, q01, q025 [, purity, ploidy]) may be NULL. Writes nothing.
// [[Rcpp::export]]
List rae_masked_means_cpp(DataFrame seg, Nullable<DataFrame> bins, Nullable<DataFrame> stats, std::string build, List params) {
  rae::Params P = params_from_list(params);
  rae::InputLog log;
  std::vector<rae::SegRow> rows = rows_from_df(seg, log);
  rae::Result R = rae::run(rows, bins_from_df(bins), stats_from_df(stats), build, P, log);
  const int nr = R.M.R, N = R.M.N;
  NumericMatrix obs((size_t)nr * N, 4), rec((size_t)nr * N, 4);
  const std::vector<char> none[2];
  std::vector<float> mean[2];
  for (int j = 0; j < N; ++j)
    for (int k0 : {rae::A0, rae::D0}) {
      rae::masked_region_means(R, (size_t)j, k0, 0, nr, none, mean);
      for (int c = 0; c < 2; ++c)
        for (int r = 0; r < nr; ++r) {
          const size_t i = (size_t)j * nr + r;
          obs(i, k0 + c) = R.M.v[k0 + c][i];
          rec(i, k0 + c) = mean[c][r];
        }
    }
  return List::create(_["observed"] = obs, _["recomputed"] = rec);
}
