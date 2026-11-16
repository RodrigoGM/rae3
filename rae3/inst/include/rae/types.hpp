// RAE 3.0 core -- shared types and parameters.
//
// Reimplementation of RAE: Taylor BS et al. (2008) PLoS ONE 3(9):e3179 (CC BY).
// Comment tags:  [STATED: where]  = the paper specifies it;
//                [INFERRED G-xx]  = interpretation, see docs/spec/02_gap_register.md.
#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace rae {

const double NaN = std::numeric_limits<double>::quiet_NaN();
enum Comp { A0 = 0, A1 = 1, D0 = 2, D1 = 3 };
static const char* const COMP_NAMES[4] = {"A0", "A1", "D0", "D1"};

// Per-sample QC / fallback flags (bit mask), reported in samples.tsv and run.json.
enum Flag : uint32_t {
  F_SD_FROM_BREAKPOINTS = 1u << 0,  // G-01a fired: per-marker SD estimated from the segmentation
  F_DN_FROM_SD          = 1u << 1,  // G-01b fired: DN derived from SD
  F_Q_FROM_PSEUDO       = 1u << 2,  // G-01c fired: D1 percentiles from pseudo-probes
  F_A0_FALLBACK         = 1u << 3,  // no gain peak: beta fixed  [STATED: S1-B]
  F_D0_FALLBACK         = 1u << 4,  // no loss peak: beta fixed  [STATED: S1-B]
  F_NO_HALFMAX          = 1u << 5,  // G-03f: half-maximum not reached inside density range
  F_D1_PCT01            = 1u << 6,  // D1 moved to the 1st percentile  [STATED: S1-B]
  F_D1_BOUND            = 1u << 7,  // D1 set to bound * E_D0          [STATED: S1-B]
  F_A1_BOUND            = 1u << 8,  // A1 set to bound * E_A0          [STATED: S1-B]
  F_UBP_EXCLUDED        = 1u << 9,  // segment-count outlier, excluded from UBP derivation
  F_SD_COHORT_MEDIAN    = 1u << 10, // G-01a: too few breakpoints, cohort median SD used
  F_STEP_FALLBACK       = 1u << 11, // G-13b: copy step unresolved; detectors placed as in RAE 2008
  F_PLOIDY_AMBIGUOUS    = 1u << 12  // G-14b: the ploidy prior's two best candidates are within ploidy_prior_ambig
};

struct Params {
  // ---- stated in the paper ------------------------------------------------------------------
  int    ubp_min_markers   = 3;      // [STATED: M&M UBP] "segments of >= 3 markers"
  int    density_n         = 2048;   // [STATED: S1-B]
  double density_from      = -1.0;   // [STATED: S1-B]
  double density_to        = 1.0;    // [STATED: S1-B]
  double density_adjust    = 2.0;    // [STATED: S1-B] "bandwidth adjustment equal to 2"
  double peak_window_lo    = 0.5;    // [STATED: S1-B] x DN ; reading [INFERRED G-03a]
  double peak_window_hi    = 2.0;    // [STATED: S1-B] x DN ; reading [INFERRED G-03a]
  double sigmoid_at_half   = 1e-3;   // [STATED: S1-B] S(H; E, beta) = 0.001
  double beta_fallback_gain = 20.0;  // [STATED: S1-B]
  double beta_fallback_loss = 20.0;  // [STATED: S1-B]; Fig. S5C shows 15 [INFERRED G-03c]
  double d1_pct_first      = 0.025;  // [STATED: S1-B]
  double d1_pct_second     = 0.01;   // [STATED: S1-B]
  double d1_bound          = 1.4;    // [STATED: S1-B]; direction [INFERRED G-04a]
  double a1_bound          = 1.4;    // [STATED: S1-B]
  double beta_a1           = 1.0 / 0.6931471805599453;  // [STATED: S1-B] 1/(log 2); natural log [INFERRED G-05a]
  int    n_perm            = 10000;  // [STATED: Main] "typically repeated 10,000 times"
  int    recleave_every    = 1000;   // [STATED: M&M]; reading [INFERRED G-08g]
  double fdr               = 0.01;   // [STATED: Main] "typical cutoff of 0.01"
  double peak_q            = 0.001;  // [STATED: M&M ROI]
  double shoulder_factor   = 2.0;    // [STATED: M&M ROI] x median analytical error
  double germline_min_score = 0.5;   // [STATED: M&M Germline]
  int    germline_min_normals = 2;   // [STATED: M&M Germline]
  double cnv_polymorphic_cover = 0.5;  // Taylor 2010: "likely polymorphic" if > 50 % of the region is
  int    cnv_polymorphic_n = 2;        //   spanned by two or more known CNVs (reported, never removed)

  // ---- inferred: every one of these is a gap-register entry -----------------------------------
  int    null_positions    = 100;    // [INFERRED G-08f] positions scored per permutation (0 = all)
  double ubp_outlier_k     = 1.5;    // [INFERRED G-07a] Tukey fence on log10(segment count); 0 disables
  bool   split_at_centromere = false;// [INFERRED G-07e]
  bool   error_sem         = false;  // [INFERRED G-10b] false: sqrt(sum/N) as rendered; true: sqrt(sum)/N
  double sd_from_breakpoints_c = 11.5;// [INFERRED G-01a] simulation-calibrated median; weak estimator, see validation report
  double marker_sd         = NaN;    // [INFERRED G-01a] explicit per-marker SD override
  double grid_spacing      = NaN;    // [INFERRED G-01d] virtual marker spacing override (bp)
  double uniform_cleave_spacing = 87000.0;  // [INFERRED G-08k] used only when no hotspot table exists
  double density_bw        = NaN;    // [INFERRED G-02a] explicit bandwidth override
  bool   density_noise_sem = true;   // [INFERRED G-02g] noise added back = SD/sqrt(n) (true) or SD (false)

  // ---- released code (RAE beta 0.91) where it differs from the clean-room reading of the text -----
  // Off under "taylor2008" and "sequencing"; preset "cgh" switches all of these. See 03_original_code_diff.md.
  bool   d0_midpoint_zero  = false;  // [CODE D-01] E_D0 = (peak + 0)/2 instead of (peak + H_d)/2
  double beta_min          = 0.0;    // [CODE D-03] |beta| of A0/D0 never below this (code: 10)
  double halfmax_min       = 0.0;    // [CODE D-04] |half-maximum point| never below this (code: 0.025)
  int    peak_span         = 3;      // [CODE D-05] local maximum over this many grid points (code: 13)
  bool   diploid_nearest_zero = false; // [CODE D-06] diploid peak = density peak nearest 0, not the global maximum
  bool   centre_window     = false;  // [CODE D-07] mode-centre on the peak of noise-free segment means within +/- DN/2
  bool   e_inflate         = false;  // [CODE D-08] grow E by 10 % until S(H; E, +/-100) <= tolerance
  bool   ubp_outlier_log   = true;   // [CODE D-09] false: boxplot rule on raw segment counts
  bool   hotspot_kde       = false;  // [CODE D-10] Gaussian KDE (bw = median gap) instead of a histogram
  bool   cleave_code       = false;  // [CODE D-11] cleavage rule as coded (left piece > half median only)
  bool   p_strict          = false;  // [CODE D-12] count null scores strictly greater than observed
  bool   roi_code          = false;  // [CODE D-13] stage-2 peaks exactly as coded
  double roi_error_factor  = 0.5;    // [CODE D-13] error is divided by this in the coded peak detector
  bool   region_weight_code = false; // [CODE D-14] within-region mean weighted by squared run length

  // ---- later-paper post-processing (NOT RAE 2008) --------------------------------------------
  double call_threshold    = 0.9;    // Taylor 2010 discrete states
  double amp_threshold     = 0.5;    // Taylor 2010 multi-copy amplification (0.25 under taylor2008 and cgh)
  double burden_threshold  = 0.75;   // rae_fga: a segment is altered when A0 or D0 >= this (Hieronymus 2014 CNA burden)
  double fga_threshold     = 0.2;    // [INFERRED G-12e] fga: altered = |log2 ratio as supplied| > this (cBioPortal 0.2)
  double fga_threshold_loss = NaN;   // [INFERRED G-12e] if set, losses use this instead (magnitude; sign ignored)
  bool   fga_assembly      = false;  // [INFERRED G-12a] rae_fga/fga denominator: false = the tumour's segment-covered span
                                     //   (first to last segment per autosome, as published); true = whole assembly

  // ---- sequencing: detectors placed from the one-copy step, not from the noise model -----------
  // On under "sequencing" (the default); off under "taylor2008" and "cgh". See 08_britroc.md.
  bool   copy_step         = true;   // [INFERRED G-13a] place E and beta from the copy-number ladder
  double comb_step_min     = 0.06;   // [INFERRED G-13b] one-copy step search range (mean-normalised ratio)
  double comb_step_max     = 0.6;    // [INFERRED G-13b]
  double comb_step_by      = 0.002;  // [INFERRED G-13b]
  double comb_rel          = 0.95;   // [INFERRED G-13b] largest step scoring >= comb_rel x the best
  double comb_tooth_min_frac = 0.02; // [INFERRED G-13c] lowest ladder level holding >= this genome fraction = 1 copy
  double amp_cn_low        = 5.0;    // [INFERRED G-13d] amplification: CN >= 5 at ploidy <= 2.7 ...
  double amp_cn_high       = 9.0;    //   ... CN >= 9 above (COSMIC / ASCAT rule, as used by BriTROC-1)
  double amp_ploidy_split  = 2.7;    // [INFERRED G-13d]
  double homdel_cn         = 0.5;    // [INFERRED G-13e] homozygous-deletion detector midway between 0 and 1 copy
  bool   a1_logistic       = true;   // [INFERRED G-13f] A1 as an ordinary sigmoid (the linear A1 is for arrays)

  // ---- opt-in analyses (off by default; default output unchanged) ------------------------------
  // Soft ploidy prior on the ladder's absolute scale: replaces the hard "lowest level = 1 copy" rule by the
  // most probable base copy number under a two-component mixture of ploidies (diploid-like, doubled).
  bool   ploidy_prior      = false;  // [INFERRED G-14a] choose the base copy number by the mixture below
  double ploidy_prior_wgd  = 0.3;    // [INFERRED G-14a] weight of the doubled component (Bielski 2018: ~30 % WGD)
  double ploidy_prior_mu1  = 1.99;   // [INFERRED G-14a] mean ploidy of non-WGD tumours (Zack 2013, TCGA ABSOLUTE)
  double ploidy_prior_mu2  = 3.31;   // [INFERRED G-14a] mean ploidy after one WGD (Zack 2013)
  double ploidy_prior_sd   = 0.4;    // [INFERRED G-14a] SD of both components; the developer's choice, NOT sourced
  double ploidy_prior_ambig = 3.0;   // [INFERRED G-14b] flag when best / runner-up density is below this

  // Per-sample baseline of each chromosome arm (the level an arm sits at when nothing is focally altered); the
  // shared base of the focal-only analysis and the arm-level test. Not used by the default analysis.
  int    arm_baseline      = 0;      // [INFERRED G-14c] 0: weighted median of the arm's segments (default); 1: two-level fit per chromosome
  int    arm_min_markers   = 50;     // [INFERRED G-14c] an arm with fewer markers is not eligible

  // Focal-only analysis: each segment is scored on its deviation from its arm's baseline, so arm-level change
  // leaves the scores, and the null is rebuilt from the same focal scores. Deviations longer than focal_max_bp are
  // scored 0 and reported per sample (focal_excluded_frac).
  bool   focal_only        = false;  // [INFERRED G-14d]
  double focal_max_bp      = 2e7;    // [INFERRED G-14d] longest focal deviation (bp); conventional upper bound 20 Mb

  // Arm-level cohort test: each sample's four detector scores at its arm baselines (G-14c), averaged over the
  // cohort per eligible arm and tested against a within-sample permutation of the arms. Writes <prefix>.arms.tsv.
  bool   arm_test          = false;  // [INFERRED G-14e]
  int    arm_n_perm        = 10000;  // [INFERRED G-14e] permutations of the arm null

  // Independence of a secondary peak: inside one event, each later stage-2 peak is re-scored with the tumours' runs
  // through the accepted peaks' apexes masked, and is independent if the residual still reaches the q <= peak_q cut.
  // Adds resid_score, resid_p and independent to roi.tsv.
  bool   peak_independence = false;  // [INFERRED G-14f]
  double peak_run_min      = 0.5;    // [INFERRED G-14f] a run is a stretch of segments scoring >= this (linear A1: amp_threshold)

  uint64_t seed            = 1;      // [INFERRED G-02b, G-11f]
  std::string preset       = "sequencing";

  // Default-constructed: "sequencing", the 2008 method with detectors placed from the copy-number
  // ladder (G-13); the discrete amplification call means "at or above the amplification copy number".
  // "taylor2008": the 2008 method as its text describes it, detectors from the noise model (array CGH).
  void apply_taylor2008() {
    preset = "taylor2008";
    copy_step = false; a1_logistic = false; amp_threshold = 0.25;
  }
  // "cgh": the released code (RAE beta 0.91), wherever it differs from the reading of the text
  // (details the text leaves open). Reproduces the published array results.
  void apply_cgh() {
    apply_taylor2008();
    preset = "cgh";
    beta_fallback_loss = 15.0;       // [CODE D-02]
    d0_midpoint_zero = true; beta_min = 10.0; halfmax_min = 0.025; peak_span = 13;
    diploid_nearest_zero = true; centre_window = true; e_inflate = true; ubp_outlier_log = false;
    hotspot_kde = true; cleave_code = true; p_strict = true; roi_code = true; region_weight_code = true;
    null_positions = 0;              // [CODE D-15] every region of every permutation is scored
    n_perm = 500; recleave_every = 50;  // [CODE D-15] PermutationCount 100 x 5 iterations, 10 cleavage sets
  }
};

// One entry of a population copy-number-variant catalogue (spec step 11).
struct CnvRow { int32_t chrom; int64_t start, end; int8_t type; float af; };  // type: +1 gain, -1 loss, 0 unknown

// Genome + ordered marker index (array probes or bins), autosomes only [INFERRED G-11a].
struct Genome {
  std::string build = "custom";
  int n_chrom = 22;
  std::vector<int64_t> chrom_len, cen_start, cen_end;  // 0 when unknown
  std::vector<std::vector<int64_t>> hotspots;          // per chromosome, sorted centres; may be empty
  bool has_hotspots = false;
  // marker grid, global index, sorted by chromosome then position
  std::vector<int32_t> m_chrom;        // 0-based chromosome index
  std::vector<int64_t> m_start, m_end; // 1-based inclusive
  std::vector<int32_t> chrom_first;    // size n_chrom + 1; markers of chromosome c are [first[c], first[c+1])
  std::vector<CnvRow> cnv;             // population CNV catalogue, sorted by chrom, start
  std::string cnv_source = "none";
  bool virtual_grid = false;           // G-01d fired
  double grid_spacing = NaN;
  int n_markers() const { return (int)m_chrom.size(); }
};

struct Segment {
  int32_t start = 0, end = 0;   // global marker indices, inclusive, after extension [INFERRED G-07d]
  int32_t chrom = 0;
  int64_t bp_start = 0, bp_end = 0;  // as given in the input
  int32_t nmark = 0;            // num.mark as given
  double mean = 0, sd = NaN;    // mean is mode-centred after normalisation
  int32_t orig = -1;            // index of the input segment it came from; set only by the centromere split  [G-14c]
};

struct Sample {
  std::string id;
  std::vector<Segment> segs;    // sorted by start
  // probe-level statistics; NaN when unknown (seg-only input)
  double dn = NaN, q01 = NaN, q025 = NaN;
  double in_purity = NaN, in_ploidy = NaN;  // optional user-supplied absolute scale  [INFERRED G-13a]
  // derived
  double sd_est = NaN;          // per-marker SD used where a segment has none
  double shift = 0;             // mode offset subtracted          [STATED: M&M]
  double bw = NaN;              // density bandwidth used
  double Ha = NaN, Hd = NaN;    // half-maximum points             [STATED: S1-B]
  double Xa = NaN, Xd = NaN;    // gain/loss peak positions used (NaN if fallback)
  double E[4] = {NaN, NaN, NaN, NaN}, B[4] = {NaN, NaN, NaN, NaN};
  uint32_t flags = 0;
  double rae_fga = NaN;         // fraction of the genome altered by the RAE scores (Hieronymus 2014 CNA burden)
  double fga = NaN;             // fraction of the genome altered by |log2 ratio| > fga_threshold
  int    n_altered = 0;         // number of segments altered by the RAE scores (as in rae_fga)
  // copy-number ladder (preset "sequencing"): source 0 none, 1 comb estimate, 2 user purity/ploidy
  int step_src = 0;
  double step = NaN, purity = NaN, ploidy = NaN, comb_score = NaN;
  double ploidy_hard = NaN, ploidy_alt = NaN, ploidy_prior_ratio = NaN;  // ploidy prior  [INFERRED G-14a]
  double focal_excluded_frac = NaN;  // focal-only: share of the covered span in deviations over focal_max_bp  [G-14d]
  double lad_m = NaN, lad_r0 = NaN;     // ladder's mean ratio and zero-copy ratio (not written to any output)
  bool a1_logistic = false;
  std::vector<double> dens_x, dens_y;  // centred density, kept for figures
};

struct Regions {                 // a breakpoint profile: regions tile every chromosome with markers
  std::vector<int32_t> start, end, chrom;  // global marker indices, inclusive
  int size() const { return (int)start.size(); }
};

inline bool is_nan(double x) { return x != x; }

}  // namespace rae
