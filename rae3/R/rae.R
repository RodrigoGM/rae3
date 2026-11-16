#' rae3: RAE reimplemented from the publication
#'
#' A reimplementation of RAE (Taylor et al. 2008, PLoS ONE 3:e3179). The algorithm lives in a
#' header-only C++ core shared with the standalone `rae` command-line tool, so both give the same
#' results for the same input and seed. Choices the paper does not determine are parameters; see
#' [rae_params()].
#'
#' @useDynLib rae3, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

#' Parameters of the RAE analysis
#'
#' Returns the full parameter list of a preset (default `"sequencing"`). Arguments override individual
#' values. Parameters marked *inferred* are not determined by the paper; each has an entry in the gap
#' register (`docs/spec/02_gap_register.md`) under the ID shown.
#'
#' @param ... named overrides.
#' @param preset `"sequencing"` (default): the RAE 2008 method with the per-tumour detectors placed on the
#'   copy-number ladder instead of from the noise model (gap register G-13): the one-copy step comes from
#'   supplied purity and ploidy (`purity_ploidy` in [rae()]) or is estimated from the periodicity of the
#'   segment means. Gain/loss detectors sit half a copy above/below the rounded ploidy, amplification at 5
#'   copies (9 above ploidy 2.7), homozygous deletion at 0.5 copies. Built for sequencing data, where the
#'   noise-model placement falls back to fixed slopes for most tumours (`docs/spec/08_britroc.md`).
#'
#'   The two array-CGH presets place the detectors from the noise model, as RAE was designed for:
#'   `"cgh"`: the released code (RAE beta 0.91), which reproduces the published array results; it fills
#'   fifteen details the text leaves open (`docs/spec/03_original_code_diff.md`). Use it for array data.
#'   `"taylor2008"`: the method as the text of Taylor et al. 2008 describes it, read without the code
#'   (`docs/spec/01_specification.md`); kept to document that reading.
#' @return A named list.
#'
#' @details
#' Stated in the paper: `ubp_min_markers` (3), `density_n` (2048), `density_from`/`density_to` (-1, 1),
#' `density_adjust` (2), `peak_window_lo`/`peak_window_hi` (0.5, 2; times DN), `sigmoid_at_half` (1e-3),
#' `beta_fallback_gain`/`beta_fallback_loss` (20), `d1_pct_first` (0.025), `d1_pct_second` (0.01),
#' `d1_bound` (1.4), `a1_bound` (1.4), `beta_a1` (1/log(2)), `n_perm` (10000), `recleave_every` (1000),
#' `fdr` (0.01), `peak_q` (0.001), `shoulder_factor` (2), `germline_min_score` (0.5),
#' `germline_min_normals` (2).
#'
#' Inferred: `null_positions` (100; G-08f), `ubp_outlier_k` (1.5; G-07a), `split_at_centromere`
#' (FALSE; G-07e), `error_sem` (FALSE; G-10b), `density_noise_sem` (TRUE; G-02g),
#' `sd_from_breakpoints_c` (11.5; G-01a), `marker_sd` (NA; G-01a), `grid_spacing` (NA; G-01d),
#' `uniform_cleave_spacing` (87000; G-08k), `density_bw` (NA; G-02a), `seed` (1).
#'
#' Released-code switches (on under `"cgh"` only): `d0_midpoint_zero`, `beta_min`, `halfmax_min`,
#' `peak_span`, `diploid_nearest_zero`, `centre_window`, `e_inflate`, `ubp_outlier_log`, `hotspot_kde`,
#' `cleave_code`, `p_strict`, `roi_code`, `roi_error_factor`, `region_weight_code`.
#'
#' Sequencing switches (on under `"sequencing"`, off under `"cgh"` and `"taylor2008"`; G-13): `copy_step`, `comb_step_min` (0.06),
#' `comb_step_max` (0.6), `comb_step_by` (0.002), `comb_rel` (0.95), `comb_tooth_min_frac` (0.02),
#' `amp_cn_low` (5), `amp_cn_high` (9), `amp_ploidy_split` (2.7), `homdel_cn` (0.5), `a1_logistic`.
#' `"sequencing"` also sets `amp_threshold` to 0.5 (0.25 under the array presets).
#'
#' From later papers, not RAE 2008: `call_threshold` (0.9) and `amp_threshold` (0.25) from
#' Taylor et al. 2010; `burden_threshold` (0.75) from Hieronymus et al. 2014.
#'
#' Fraction of the genome altered per sample (see [rae()], section *Output*): `burden_threshold` (0.75)
#' for `rae_fga`, `fga_threshold` (0.2) and `fga_threshold_loss` (`NULL`: the same) for `fga`, and `fga_assembly` (`FALSE`) for the denominator of
#' both (gap register G-12a, G-12e). The same under every preset.
#'
#' Opt-in analyses (off by default, so default output is unchanged): `ploidy_prior` (`FALSE`; G-14a)
#' replaces the rule that the lowest ladder level holding `comb_tooth_min_frac` of the genome is one copy
#' by a soft prior. Each admissible base copy number, from that rule's choice up to the largest purity <= 1
#' allows, gives a mean ploidy; the one with the highest density under a mixture of two normal
#' distributions is taken. The mixture has weight `ploidy_prior_wgd` (0.3; G-14a; about 30 % of advanced
#' cancers are genome-doubled, Bielski et al. 2018 Nat Genet 50:1189) on a doubled-genome component with mean
#' `ploidy_prior_mu2` (3.31; G-14a), the rest on a diploid-like component with mean `ploidy_prior_mu1`
#' (1.99; G-14a; the two means are the mean ploidy of non-doubled tumours and of tumours after one doubling in
#' Zack et al. 2013 Nat Genet 45:1134), both with standard deviation `ploidy_prior_sd` (0.4; G-14a; the
#' developer's choice, not sourced). `ploidy_prior_mu1` and `ploidy_prior_mu2` must be finite, `ploidy_prior_sd` > 0,
#' `ploidy_prior_wgd` between 0 and 1 and `ploidy_prior_ambig` >= 1, or the run stops. `ploidy_prior_ratio` is `NA` where it would not be finite (the
#' runner-up's density underflows). If even the best candidate's density underflows to 0 the prior gave no
#' information: the hard choice is kept, the ratio is `NA` and `ploidy_ambiguous` is set. It can only
#' raise the base copy number, and applies only when the ladder is estimated (not with supplied purity and
#' ploidy). `ploidy_prior_ambig` (3; G-14b) sets the flag `ploidy_ambiguous`, raised when the best candidate
#' is less than this many times as probable as the runner-up. With the prior on, `samples` gains the columns
#' `ploidy_hard`, `ploidy_alt` and `ploidy_prior_ratio`; these are `NA` for tumours where the prior did not
#' run (a single admissible base copy number, supplied purity and ploidy, or a fallback ladder). On BriTROC-1 the
#' pre-registered benchmark found no improvement from the prior (`docs/spec/08_britroc.md`, section 8).
#'
#' Arm baselines (G-14c; the shared base of later opt-in analyses, so they change no default output):
#' `arm_baseline` (`"median"` or `"twolevel"`) says how the level an arm sits at is estimated. `"median"` (the
#' default, chosen by a pre-registered comparison on simulated data and BriTROC-1) is the marker-weighted median of the arm's segments (the lower one: the smallest value at which the cumulative
#' weight reaches half). `"twolevel"` fits each chromosome one level or two, split at a single breakpoint, by
#' least absolute deviations and the Bayesian information criterion with penalty (2k - 1) ln n for k levels (the
#' broad-level fit described in the Supplementary Methods of Mermel et al. 2011 Genome Biol 12:R41), and each arm takes the level
#' that covers most of its markers. `arm_min_markers` (50) is the number of markers an arm needs to be eligible.
#'
#' Focal-only analysis (G-14d): `focal_only` (`FALSE`) scores each segment on its deviation from its arm's baseline,
#' so arm-level change leaves the scores and the significance test, and the null model is rebuilt from the same
#' focal scores. Segments crossing the centromere are first cut there (a segment with fewer than 2 markers is left
#' whole and judged on the arm of its start marker). The gain, amplification and loss detectors are evaluated at the residual: in copy units
#' kmode + (c - b) copies, kmode being the sample's rounded ploidy (at least 1), c the segment's and b the arm's copy
#' number; on the log2 scale (arrays, and samples whose copy step fell back) the segment mean minus b. The
#' homozygous-deletion detector stays on the absolute value, so a homozygous deletion on a lost arm keeps it. A run of
#' contiguous segments deviating in one direction over more than `focal_max_bp` (2e7, 20 Mb, the conventional upper
#' bound of a focal event) is scored 0 and its share of the sample's segment-covered span is written as
#' `focal_excluded_frac` in `samples`; such a run continues across the centromere, so a deviation spanning both arms
#' is one event. Focal (split segments, focal scores): the breakpoint profile, the null model, `regions`, `roi`,
#' `calls`, `segments` and `n_segments`. Absolute, and equal to a run without `focal_only`: `rae_fga`, `fga` and `n_altered` (the cut keeps the
#' segment's bp length, and a cut pair is counted once). `focal_max_bp` must be > 0. With `build =
#' "custom"` (no centromeres) each chromosome is one arm, and the run records a warning.
#'
#' Arm-level cohort test (G-14e): `arm_test` (`FALSE`) tests whether each chromosome arm is altered in more
#' tumours than a random arm, given each tumour's own number of altered arms. Each tumour's four detectors are
#' evaluated at its arm baselines (in copy units through the ladder, otherwise on the log2 scale), the scores are
#' averaged over the tumours with data on the arm and combined as for regions (gain sqrt(A0^2 + A1^2), loss
#' sqrt(D0^2 + D1^2)). The null permutes each tumour's arm scores across the eligible arms (at least
#' `arm_min_markers` markers) `arm_n_perm` (10000) times, on its own random stream. The table, written as a ninth
#' file `<prefix>.arms.tsv` and returned as `arms`, has pooled permutation p-values, Benjamini-Hochberg q and
#' Westfall-Young max-T adjusted p-values (`fwer_*`) for gain and loss separately. In highly altered cohorts
#' (high-grade serous ovarian cancer, for example) the null mean is high, so the test finds the arms altered
#' more than the tumours' average arm. It changes no other output. `arm_n_perm` and `arm_min_markers` must be
#' at least 1. With `build = "custom"` each chromosome is one arm, and the run records a warning.
#' On the command line `--focal-max-mb` is in Mb while `focal_max_bp` is in bp, and `--arm-perm` sets `arm_n_perm`.
#'
#' Independence of a secondary peak (G-14f): `peak_independence` (`FALSE`) asks, for an event with more than one
#' stage-2 peak, whether each later peak is more than the shoulder of a higher one. The peaks are taken in order of
#' score; the highest is the primary. For each later peak, every tumour's runs (maximal stretches of consecutive
#' segments on one chromosome whose detector score is at least `peak_run_min`, 0.5; for the linear A1 of the array
#' presets and of samples flagged `step_fallback`, `amp_threshold`) that reach the apex region of the primary or of a peak already found independent are set
#' to 0, component by component (A0 and A1 for gain, D0 and D1 for loss). The tumours' region means over the peak's
#' span are recomputed with the same weighting, and the residual summary score is their cohort mean combined as for
#' regions; `resid_score` is its maximum over the span and `resid_p` its p-value against the same null as the
#' observed scores. The peak is `independent` (1) when `resid_score` is at least the lowest observed score among
#' regions with q <= `peak_q`, and then joins the accepted set; otherwise 0. The primary peak, a peak alone in its
#' event and the event rows are `NA`. Under `focal_only` the runs use the focal scores, so they describe focal
#' events only, and they continue across the centromere. The idea is GISTIC2's peel-off (Mermel et al. 2011 Genome
#' Biol 12:R41), reimplemented from its published description; no GISTIC2 code is used. With the switch on, `roi`
#' gains the columns `resid_score`, `resid_p` and `independent`; nothing else changes. `peak_run_min` must be in
#' (0, 1]. Sub-threshold stretches (score below `peak_run_min`) are never removed, so a shoulder carried by weak
#' runs can still be called independent; raising `peak_run_min` is the control.
#' @export
rae_params <- function(..., preset = c("sequencing", "cgh", "taylor2008")) {
  preset <- match.arg(preset)
  p <- list(
    ubp_min_markers = 3L, density_n = 2048L, density_from = -1, density_to = 1, density_adjust = 2,
    peak_window_lo = 0.5, peak_window_hi = 2, sigmoid_at_half = 1e-3,
    beta_fallback_gain = 20, beta_fallback_loss = 20, d1_pct_first = 0.025, d1_pct_second = 0.01,
    d1_bound = 1.4, a1_bound = 1.4, beta_a1 = 1 / log(2), n_perm = 10000L, recleave_every = 1000L,
    fdr = 0.01, peak_q = 0.001, shoulder_factor = 2, germline_min_score = 0.5, germline_min_normals = 2L,
    null_positions = 100L, ubp_outlier_k = 1.5, split_at_centromere = FALSE, error_sem = FALSE,
    density_noise_sem = TRUE, sd_from_breakpoints_c = 11.5, marker_sd = NULL, grid_spacing = NULL,
    uniform_cleave_spacing = 87000, density_bw = NULL,
    d0_midpoint_zero = FALSE, beta_min = 0, halfmax_min = 0, peak_span = 3L, diploid_nearest_zero = FALSE,
    centre_window = FALSE, e_inflate = FALSE, ubp_outlier_log = TRUE, hotspot_kde = FALSE,
    cleave_code = FALSE, p_strict = FALSE, roi_code = FALSE, roi_error_factor = 0.5,
    region_weight_code = FALSE,
    copy_step = FALSE, comb_step_min = 0.06, comb_step_max = 0.6, comb_step_by = 0.002, comb_rel = 0.95,
    comb_tooth_min_frac = 0.02, amp_cn_low = 5, amp_cn_high = 9, amp_ploidy_split = 2.7, homdel_cn = 0.5,
    a1_logistic = FALSE,
    ploidy_prior = FALSE, ploidy_prior_wgd = 0.3, ploidy_prior_mu1 = 1.99, ploidy_prior_mu2 = 3.31,
    ploidy_prior_sd = 0.4, ploidy_prior_ambig = 3,
    arm_baseline = "median", arm_min_markers = 50L, focal_only = FALSE, focal_max_bp = 2e7,
    arm_test = FALSE, arm_n_perm = 10000L, peak_independence = FALSE, peak_run_min = 0.5,
    call_threshold = 0.9, amp_threshold = 0.25, burden_threshold = 0.75, fga_threshold = 0.2,
    fga_threshold_loss = NULL, fga_assembly = FALSE, seed = 1)
  if (preset == "cgh") {
    orig <- list(beta_fallback_loss = 15, d0_midpoint_zero = TRUE, beta_min = 10, halfmax_min = 0.025,
                 peak_span = 13L, diploid_nearest_zero = TRUE, centre_window = TRUE, e_inflate = TRUE,
                 ubp_outlier_log = FALSE, hotspot_kde = TRUE, cleave_code = TRUE, p_strict = TRUE,
                 roi_code = TRUE, region_weight_code = TRUE, null_positions = 0L, n_perm = 500L,
                 recleave_every = 50L)
    p[names(orig)] <- orig
  }
  if (preset == "sequencing") p[c("copy_step", "a1_logistic", "amp_threshold")] <- list(TRUE, TRUE, 0.5)
  o <- list(...)
  bad <- setdiff(names(o), names(p))
  if (length(bad)) stop("unknown RAE parameter(s): ", paste(bad, collapse = ", "))
  for (n in names(o)) p[n] <- list(o[[n]])
  # the core takes the method as a code: 1 two-level fit, 0 weighted median
  if (is.character(p$arm_baseline)) p$arm_baseline <- c(median = 0L, twolevel = 1L)[[match.arg(p$arm_baseline, c("median", "twolevel"))]]
  if (!isTRUE(p$arm_baseline %in% 0:1)) stop("arm_baseline must be \"twolevel\" or \"median\"")
  p$preset <- preset
  p
}

.seg_df <- function(seg) {
  seg <- as.data.frame(seg, stringsAsFactors = FALSE)
  if (ncol(seg) < 6) stop("a seg table needs six columns: ID, chrom, loc.start, loc.end, num.mark, seg.mean")
  out <- data.frame(ID = as.character(seg[[1]]), chrom = as.character(seg[[2]]),
                    loc.start = as.numeric(seg[[3]]), loc.end = as.numeric(seg[[4]]),
                    num.mark = as.numeric(seg[[5]]), seg.mean = as.numeric(seg[[6]]), stringsAsFactors = FALSE)
  if (ncol(seg) >= 7) out$seg.sd <- as.numeric(seg[[7]])
  out
}
.bins_df <- function(bins) {
  if (is.null(bins)) return(NULL)
  bins <- as.data.frame(bins, stringsAsFactors = FALSE)
  data.frame(chrom = as.character(bins[[1]]), start = as.numeric(bins[[2]]), end = as.numeric(bins[[3]]),
             stringsAsFactors = FALSE)
}
.stats_df <- function(s, purity_ploidy = NULL) {
  if (!is.null(purity_ploidy)) {
    pp <- as.data.frame(purity_ploidy, stringsAsFactors = FALSE)
    if (!all(c("ID", "purity", "ploidy") %in% names(pp))) stop("purity_ploidy needs columns: ID, purity, ploidy")
    pp <- data.frame(ID = as.character(pp$ID), purity = as.numeric(pp$purity), ploidy = as.numeric(pp$ploidy),
                     stringsAsFactors = FALSE)
    s <- if (is.null(s)) pp else merge(as.data.frame(s, stringsAsFactors = FALSE)[, setdiff(names(s), c("purity", "ploidy")), drop = FALSE],
                                       pp, by = "ID", all = TRUE, sort = FALSE)
  }
  if (is.null(s)) return(NULL)
  s <- as.data.frame(s, stringsAsFactors = FALSE)
  if (!"ID" %in% names(s)) stop("sample statistics need an ID column")
  has_noise <- any(c("dn", "q01", "q025") %in% names(s))
  if (has_noise && !all(c("dn", "q01", "q025") %in% names(s))) stop("sample statistics need columns: ID, dn, q01, q025")
  num <- function(n) if (n %in% names(s)) as.numeric(s[[n]]) else rep(NA_real_, nrow(s))
  out <- data.frame(ID = as.character(s$ID), dn = num("dn"), q01 = num("q01"), q025 = num("q025"), stringsAsFactors = FALSE)
  if (any(c("purity", "ploidy") %in% names(s))) { out$purity <- num("purity"); out$ploidy <- num("ploidy") }
  out
}

#' Run RAE
#'
#' Maps recurrent copy-number alterations in a cohort: per-tumour detectors, unified breakpoint
#' profile, permutation null model, q-values and regions of interest.
#'
#' @param x A `QDNAseqCopyNumbers` object that has been segmented (`QDNAseq::segmentBins`), or a
#'   seg `data.frame` with columns ID, chrom, loc.start, loc.end, num.mark, seg.mean (log2 ratio)
#'   and optionally seg.sd, in that order.
#' @param ... passed to methods.
#' @return An object of class `rae_result`: a list with data frames `regions`, `roi`, `samples`,
#'   `segments`, `calls`, the `manifest` (parsed `run.json` text), `params`, and `files` (the default
#'   output set, identical to what the command-line tool writes); with `arm_test` in [rae_params()] also `arms`,
#'   the arm-level table. See *Output*.
#'
#' @section Output:
#' `samples` has one row per tumour. Besides the fitted detectors it reports how much of the genome is
#' altered, as fractions (0.054 = 5.4 %). CNA burden is the concept; the fraction of the genome altered
#' (FGA) is how it is measured. rae3 gives two FGAs, which differ in what counts as altered:
#' \describe{
#'   \item{`rae_fga`}{FGA by the RAE calls: the territory of segments whose gain or loss score A0 or D0
#'     is at least `burden_threshold` (0.75). This is the CNA burden of Taylor et al. 2010 and
#'     Hieronymus et al. 2014, whose 5.41 % cut-off separates prostate cancers that recur more often. It
#'     follows RAE's own assessment of each tumour (its noise, and under `"sequencing"` its copy-number
#'     ladder), so a gain or loss smaller than any fixed log2 cut-off still counts, and it depends on the
#'     preset.}
#'   \item{`fga`}{FGA by a fixed threshold: the territory of segments whose log2 ratio, as supplied in
#'     the seg input, is above `fga_threshold` (0.2) or below minus `fga_threshold_loss` (by default the
#'     same 0.2), whatever their RAE scores. With 0.2 this
#'     is the FGA of cBioPortal and MSK-IMPACT ("log2 copy number gain or loss > 0.2 divided by the size
#'     of the genome profiled"). 0.2 is a convention for bulk tumours, where normal cells dilute each
#'     copy-number step towards 0. Without that dilution (single-cell data, purity 1, diploid) one copy
#'     gained is log2(3/2) = 0.58 and one copy lost is log2(1/2) = -1, so the thresholds are set higher,
#'     each just short of its one-copy step with a buffer for noise. Because the two steps differ,
#'     `fga_threshold_loss` sets the loss side separately (as a magnitude: 0.8 and -0.8 both mean
#'     log2 ratio < -0.8); left `NULL`, losses use `fga_threshold`. It uses neither the RAE scores nor
#'     RAE's re-centring of each tumour, so it is the same under every preset.}
#'   \item{`n_altered`}{the number of segments counted in `rae_fga`, not divided by anything. With `focal_only`
#'     it counts input segments: a segment cut at the centromere counts once if either piece is altered.}
#' }
#' Both fractions use the same denominator, set by `fga_assembly` in [rae_params()] (gap register G-12a):
#' \describe{
#'   \item{`FALSE` (default, every preset): the segment-covered span}{the tumour's measured genome: on each
#'     autosome, from the start of its first segment to the end of its last, summed. This is the published
#'     definition: `rae_fga` reproduces Hieronymus et al. 2014 and Taylor et al. 2010 to a ratio of 0.998,
#'     so the 5.41 % cut-off applies as published. It leaves out what no assay measured (chromosome ends,
#'     the short arms of chromosomes 13, 14, 15, 21 and 22: about 3 % of the assembly on arrays) and
#'     follows each platform's coverage, so it suits array CGH and sequencing alike.}
#'   \item{`TRUE`: the whole assembly}{the summed autosome lengths of the reference build, the same for
#'     every tumour. A fixed yardstick across assays with different coverage, but lower by the unmeasured
#'     share (about 0.967 x the published values on arrays).}
#' }
#' Output of the opt-in analyses (see [rae_params()]; each is off by default, and default output is unchanged):
#' \describe{
#'   \item{`arms.tsv` (`arm_test`)}{a ninth file, `<prefix>.arms.tsv`, written only when the arm-level test
#'     is on (any older file under the prefix is removed otherwise). One row per eligible arm: `arm`,
#'     `chrom`, `arm_label` (`1p`, `1q`; a bare chromosome number such as `1` where the build has no
#'     centromere, so the whole chromosome is one arm), `start`, `end`, `length_bp`, `n_markers`,
#'     `n_samples` (tumours with data on the arm), the mean detector scores `A0 A1 D0 D1`, `score_gain` and
#'     `score_loss`, and for gain and loss in turn the permutation p-value (`p_*`), Benjamini-Hochberg q
#'     (`q_*`) and Westfall-Young max-T adjusted p (`fwer_*`). `frac_gain` and `frac_loss` are the shares of
#'     tumours with a gain or loss call at `call_threshold`; `baseline_coverage` is the mean share of the
#'     arm's markers at the tumour's arm baseline. `run.json` records `arms_significant_gain` and
#'     `arms_significant_loss` (q <= `fdr`).}
#'   \item{`roi.tsv` (`peak_independence`)}{three columns appended: `resid_score` (the peak's score after
#'     removing the alterations that carry the higher peaks of the event), `resid_p` (its empirical p against
#'     the same null; a heuristic comparison with the observed cut, not a calibrated p-value) and
#'     `independent` (1, 0 or `NA`). `NA` for the primary peak, a peak alone in its event, and stage-1 rows.}
#'   \item{`samples.tsv` (`ploidy_prior`)}{`ploidy_hard` (mean ploidy by the lowest-level rule),
#'     `ploidy_alt` (the runner-up candidate) and `ploidy_prior_ratio` (density of the chosen candidate over
#'     the runner-up's). All three are `NA` for tumours with a single admissible base copy number, with
#'     supplied purity and ploidy, and for fallback ladders. The `flags` column then can hold
#'     `ploidy_ambiguous`, set when the ratio is below `ploidy_prior_ambig` (3), or when the prior's densities
#'     underflowed and the hard choice was kept.}
#'   \item{`samples.tsv` (`focal_only`)}{`focal_excluded_frac`: the share of the tumour's segment-covered span
#'     in runs longer than `focal_max_bp`, whose scores were set to 0. `regions`, `roi`, `calls` and
#'     `segments` then describe focal events; `rae_fga`, `fga` and `n_altered` stay absolute.}
#' }
#' Before version 3.0.0.9000 the RAE-based fraction was named `cna_burden` and divided by the whole
#' assembly; [rae_read()] renames it `rae_fga` and sets `fga_assembly = TRUE` for runs saved then (they
#' have no threshold-based `fga`).
#' @export
rae <- function(x, ...) UseMethod("rae")

#' @rdname rae
#' @param bins optional marker table (chrom, start, end). Without it a uniform marker grid is
#'   reconstructed from the segments (gap G-01d).
#' @param sample_stats optional data frame with columns ID, dn, q01, q025 computed from probe-level
#'   log2 ratios (derivative noise; 1st and 2.5th percentiles). Without it they are derived from the
#'   segments (gaps G-01b, G-01c).
#' @param normals,normal_stats optional seg table and statistics of non-tumour samples, used only to
#'   flag germline regions.
#' @param purity_ploidy optional data frame with columns ID, purity (0-1] and ploidy (mean copy number),
#'   e.g. from ACE, Rascal or ichorCNA. Used only by `preset = "sequencing"` to place the detectors on the
#'   copy-number ladder; tumours without a row are estimated from the data (G-13). Can also be given as
#'   `purity`/`ploidy` columns of `sample_stats`.
#' @param build `"hg38"` (default), `"hg19"`, `"chm13"`, `"auto"` or `"custom"`. (`"hg17"`, the
#'   paper's build, is still accepted for development work on paper-era data and raises a warning
#'   in the manifest.)
#' @param params list from [rae_params()].
#' @param cnv population copy-number-variant catalogue used for the overlap columns of `roi`:
#'   `"default"` (bundled: DGV Gold Standard for hg19, gnomAD-SV v4.1 for hg38; HPRC long-read deletions plus lifted gnomAD gains for
#'   CHM13), `"none"`, or the path of a table with columns chrom, start, end, type (`gain`/`loss`)
#'   and optionally allele frequency. Overlap is reported; nothing is filtered.
#' @param out output prefix. Defaults to a temporary location.
#' @param label text recorded as the input name in the manifest and report.
#' @param filter_germline_overlaps what to do with events that overlap germline variation, defined as
#'   `likely_polymorphic` (more than half covered by at least two population CNVs) or seen in at least
#'   `germline_min_normals` of the supplied normals. `FALSE` (default, as in the paper: reported, never
#'   removed): kept in `roi`, flagged in the column `germline_overlap` and drawn grey with the mark "Г" by the
#'   plotting functions. `TRUE`: moved from `roi` to `roi_germline`. q-values are never changed. The files
#'   written to `out` are always the full, unfiltered set.
#' @export
rae.data.frame <- function(x, bins = NULL, sample_stats = NULL, normals = NULL, normal_stats = NULL,
                           purity_ploidy = NULL, build = "hg38", params = rae_params(), out = NULL, label = "data.frame",
                           cnv = "default", filter_germline_overlaps = FALSE, ...) {
  if (is.null(out)) out <- file.path(tempfile("rae_"), "cohort.rae")
  run <- rae_run_cpp(.seg_df(x), .bins_df(bins), .stats_df(sample_stats, purity_ploidy),
                     if (is.null(normals)) NULL else .seg_df(normals), .stats_df(normal_stats),
                     build, params, path.expand(out), label, cnv)
  # the result is what the files say, exactly as for a command-line run; only the parameter list is the
  # caller's own object (same values, R types)
  res <- rae_read(path.expand(out), filter_germline_overlaps = filter_germline_overlaps)
  res$files <- run$files
  res$params <- params
  res
}

#' Read the output of the command-line tool
#'
#' Loads the default output set written by `rae` (the command-line tool) or by [rae()] into the object [rae()]
#' returns, so every plotting function works on command-line results. `rae()` itself is built on this reader,
#' so both give the same object for the same files.
#'
#' @param prefix output prefix, as given to `rae --out` (e.g. `"results/cohort.rae"`).
#' @param filter_germline_overlaps as in [rae()].
#' @return a `rae_result`. Parameters come from `run.json`, typed as in [rae_params()]. The null-model
#'   internals used by [rae_plot_hotspots()] and [rae_plot_region_sizes()] come from `<prefix>.null.tsv`;
#'   outputs written before that file existed are read without them (with a message). With `arm_test`, the
#'   arm-level table `<prefix>.arms.tsv` is read into `arms` (used by [rae_plot_arms()]; `arm_label` is a bare
#'   chromosome number where the build has no centromere); otherwise (`arm_test` off in `run.json`) `arms` is `NULL`. With `peak_independence`, the `roi` columns `resid_score` and `resid_p` are
#'   numeric and `independent` is integer (0, 1, `NA`). The `samples` columns `ploidy_hard`, `ploidy_alt`,
#'   `ploidy_prior_ratio` (`ploidy_prior`) and `focal_excluded_frac` (`focal_only`) are numeric.
#' @export
rae_read <- function(prefix, filter_germline_overlaps = FALSE) {
  f <- function(suffix) paste0(prefix, suffix)
  need <- f(c(".regions.tsv", ".roi.tsv", ".samples.tsv", ".segments.tsv", ".calls.tsv", ".run.json"))
  miss <- need[!file.exists(need)]
  if (length(miss)) stop("RAE output not found: ", paste(miss, collapse = ", "), call. = FALSE)
  rd <- function(suffix) utils::read.delim(f(suffix), check.names = FALSE, stringsAsFactors = FALSE)
  js <- readLines(f(".run.json"), warn = FALSE)
  man <- .read_manifest(js)
  diag <- list(build = man$build)
  if (file.exists(f(".null.tsv"))) diag <- .read_null(f(".null.tsv"), man)
  else message("no ", basename(f(".null.tsv")), " (output of an older version): rae_plot_hotspots() and ",
               "rae_plot_region_sizes() are unavailable for this result")
  files <- f(c(".regions.tsv", ".roi.tsv", ".samples.tsv", ".segments.tsv", ".calls.tsv", ".null.tsv", ".run.json", ".report.html",
               if (isTRUE(man$params$arm_test)) ".arms.tsv"))
  smp <- rd(".samples.tsv")
  roi <- rd(".roi.tsv")
  # independence of a secondary peak (peak_independence only, G-14f): typed even when every value is NA
  for (n in intersect(c("resid_score", "resid_p"), names(roi))) roi[[n]] <- as.numeric(roi[[n]])
  # ploidy prior (G-14a, G-14b) and focal-only (G-14d) columns: numeric even when every value is NA or 0
  for (n in intersect(c("ploidy_hard", "ploidy_alt", "ploidy_prior_ratio", "focal_excluded_frac"), names(smp)))
    smp[[n]] <- as.numeric(smp[[n]])
  if (!is.null(roi$independent)) roi$independent <- as.integer(roi$independent)
  names(smp)[names(smp) == "cna_burden"] <- "rae_fga"   # name before 3.0.0.9000 (assembly denominator)
  res <- structure(list(regions = rd(".regions.tsv"), roi = roi, samples = smp,
                        segments = rd(".segments.tsv"), calls = rd(".calls.tsv"), manifest = js,
                        params = man$params, files = files[file.exists(files)], build = man$build, diagnostics = diag,
                        measured_noise = man$measured_noise),
                   class = "rae_result")
  # the arm-level test (arm_test only, G-14e)
  # only when this run had the test on, so a file left by an earlier run under the same prefix is never read
  if (isTRUE(man$params$arm_test) && file.exists(f(".arms.tsv"))) res$arms <- utils::read.delim(f(".arms.tsv"), check.names = FALSE, stringsAsFactors = FALSE,
                                                                 colClasses = c(arm_label = "character"))
  rae_germline(res, filter = filter_germline_overlaps)
}

# run.json as written by report.hpp: the build, the hotspot flag, whether noise was measured, and every
# parameter, typed after rae_params() of the same preset (JSON null -> NULL, 0/1 -> logical where logical).
.read_manifest <- function(js) {
  txt <- paste(js, collapse = "\n")
  obj <- function(name) {
    m <- regmatches(txt, regexpr(paste0('"', name, '": \\{[^{}]*\\}'), txt))
    if (!length(m)) stop("run.json has no '", name, "' block", call. = FALSE)
    m <- sub('^"[^"]+": \\{', "", m)   # the block's contents only
    kv <- regmatches(m, gregexpr('"[^"]+": ("[^"]*"|[^,}[:space:]]+)', m))[[1]]
    stats::setNames(lapply(sub('^"[^"]+": ', "", kv), function(v) if (grepl('^"', v)) gsub('"', "", v) else v),
                    sub('^"([^"]+)".*', "\\1", kv))
  }
  g <- obj("genome"); pj <- obj("parameters"); fb <- obj("inferred_fallbacks_fired")
  legacy <- c(paper = "taylor2008", original = "cgh")       # preset names before rae3 3.0.0.9000
  if (pj$preset %in% names(legacy)) pj$preset <- legacy[[pj$preset]]
  p <- rae_params(preset = pj$preset)
  if (is.null(pj$fga_assembly)) pj$fga_assembly <- "1"   # runs before the switch divided by the assembly
  for (n in setdiff(intersect(names(pj), names(p)), "preset")) {   # core-only constants are not R parameters
    v <- pj[[n]]
    p[n] <- list(if (identical(v, "null")) NULL
                 else if (is.logical(p[[n]])) as.numeric(v) != 0
                 else if (is.integer(p[[n]])) as.integer(v)
                 else as.numeric(v))
  }
  list(build = g$build, has_hotspots = identical(g$hotspot_table, "true"), params = p,
       measured_noise = identical(fb$`G-01b_dn_from_sd_samples`, "0"))
}

# <prefix>.null.tsv (item, chrom, value) -> the diagnostics list of a rae_result.
.read_null <- function(path, man) {
  nl <- utils::read.delim(path, colClasses = c("character", "integer", "numeric"))
  pick <- function(item) nl[nl$item == item, , drop = FALSE]
  per_chrom <- function(item) { z <- pick(item); z$value[order(z$chrom)] }
  pts <- function(item) { z <- pick(item); data.frame(chrom = z$chrom, pos = z$value) }
  list(hotspots = pts("hotspot"), sampled_cuts = pts("sampled_cut"),
       chrom_len = per_chrom("chrom_len"), cen_start = per_chrom("cen_start"), cen_end = per_chrom("cen_end"),
       ubp_bp = pick("ubp_bp")$value, cleaved_bp = pick("cleaved_bp")$value,
       has_hotspots = man$has_hotspots, build = man$build)
}

#' Flag or remove events that overlap germline variation
#'
#' Adds the logical column `germline_overlap` to `x$roi`: `TRUE` when the event is `likely_polymorphic`
#' (more than half of it covered by at least two population CNVs) or was seen in at least
#' `germline_min_normals` supplied normal samples. With `filter = TRUE` those rows are moved to
#' `x$roi_germline`; with `FALSE` they stay and the plotting functions draw them grey with the mark "Г".
#' q-values are not recomputed.
#'
#' @param x a `rae_result`.
#' @param filter logical, see above.
#' @return the modified `rae_result`.
#' @export
rae_germline <- function(x, filter = FALSE) {
  if (!is.null(x$roi_germline)) { x$roi <- rbind(x$roi, x$roi_germline); x$roi <- x$roi[order(x$roi$roi), , drop = FALSE] }
  r <- x$roi
  poly <- if (is.null(r$likely_polymorphic)) FALSE else r$likely_polymorphic %in% c(1, TRUE)
  nrm <- if (is.null(r$germline_normals_max)) FALSE else r$germline_normals_max >= x$params$germline_min_normals
  r$germline_overlap <- poly | nrm
  x$roi <- r[!(filter & r$germline_overlap), , drop = FALSE]
  x$roi_germline <- if (filter) r[r$germline_overlap, , drop = FALSE] else NULL
  x$filter_germline_overlaps <- isTRUE(filter)
  x
}

#' @rdname rae
#' @param bin_noise `FALSE` (default): noise is estimated from the segmentation, exactly as the
#'   command-line tool does for a plain `.seg`, so both give the same result. `TRUE`: the bin-to-bin
#'   noise (the paper's DN, the same statistic as MAPD), the 1st and 2.5th percentiles and the segment
#'   standard deviations are measured from the unsegmented `copynumber` values the object already holds.
#' @export
rae.QDNAseqCopyNumbers <- function(x, normals = NULL, build = "hg38", params = rae_params(), out = NULL,
                                   label = "QDNAseqCopyNumbers", bin_noise = FALSE, ...) {
  inp <- as_seg(x, bin_noise = bin_noise)
  nrm <- if (is.null(normals)) NULL else if (is.data.frame(normals)) list(seg = normals, stats = NULL) else
    as_seg(normals, bin_noise = bin_noise)
  rae.data.frame(inp$seg, bins = inp$bins, sample_stats = inp$stats, normals = nrm$seg, normal_stats = nrm$stats,
                 build = build, params = params, out = out, label = label, ...)
}

#' Convert a segmented QDNAseq object to RAE input
#'
#' Run-length encodes the `segmented` assay into a seg table. With `bin_noise = TRUE` it also computes
#' from the bin-level `copynumber` assay the statistics the paper derives from probe-level data: the
#' segment standard deviation, the derivative noise DN (Eq. 1 of the paper; MAPD) and the 1st and
#' 2.5th percentiles.
#' Only autosomal bins with `use == TRUE` are kept; a bin that is `NA` for a sample is skipped for
#' that sample.
#'
#' @param x a `QDNAseqCopyNumbers` object with `copynumber` and `segmented` assay data.
#' @param prefix if given, also writes `<prefix>.seg`, `<prefix>.bins.tsv` and, with `bin_noise = TRUE`,
#'   `<prefix>.sample_stats.tsv` for the command-line tool.
#' @param bin_noise `FALSE` (default): plain six-column seg and no sample statistics. `TRUE`: adds
#'   `seg.sd` and the per-sample statistics measured from the unsegmented bins.
#' @param transform `"log2"` (QDNAseq stores ratios) or `"none"` if the assays already hold log2 ratios.
#' @return list with `seg`, `bins`, `stats` (`NULL` unless `bin_noise = TRUE`).
#' @export
as_seg <- function(x, prefix = NULL, transform = c("log2", "none"), bin_noise = FALSE) {
  transform <- match.arg(transform)
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("package 'Biobase' is needed for QDNAseq objects")
  fd <- Biobase::fData(x)
  cn <- Biobase::assayDataElement(x, "copynumber")
  if (!"segmented" %in% Biobase::assayDataElementNames(x)) stop("object is not segmented; run QDNAseq::segmentBins() first")
  sg <- Biobase::assayDataElement(x, "segmented")
  chrom <- sub("^chr", "", as.character(fd$chromosome))
  keep <- chrom %in% as.character(1:22) & (if (is.null(fd$use)) TRUE else fd$use)
  ord <- order(as.integer(chrom[keep]), fd$start[keep])
  idx <- which(keep)[ord]
  chrom <- as.integer(chrom[idx]); start <- fd$start[idx]; end <- fd$end[idx]
  tf <- if (transform == "log2") function(m) log2(pmax(m, 0) + .Machine$double.xmin) else identity
  cn <- tf(cn[idx, , drop = FALSE]); sg <- tf(sg[idx, , drop = FALSE])
  ids <- colnames(cn)
  segs <- vector("list", length(ids)); stats <- vector("list", length(ids))
  for (j in seq_along(ids)) {
    ok <- which(!is.na(sg[, j]) & !is.na(cn[, j]) & is.finite(cn[, j]))
    v <- sg[ok, j]; ch <- chrom[ok]; y <- cn[ok, j]
    brk <- c(TRUE, v[-1] != v[-length(v)] | ch[-1] != ch[-length(ch)])
    run <- cumsum(brk)
    first <- which(brk); last <- c(first[-1] - 1L, length(v))
    n <- last - first + 1L
    sdv <- as.numeric(tapply(y, run, stats::sd)); sdv[is.na(sdv)] <- 0
    segs[[j]] <- data.frame(ID = ids[j], chrom = ch[first], loc.start = start[ok][first], loc.end = end[ok][last],
                            num.mark = n, seg.mean = v[first], seg.sd = sdv, stringsAsFactors = FALSE)
    same <- ch[-1] == ch[-length(ch)]                       # DN within chromosomes  [gap G-11b]
    q <- stats::quantile(y, c(0.01, 0.025), names = FALSE)   # R type 7              [gap G-04c]
    stats[[j]] <- data.frame(ID = ids[j], dn = stats::median(abs(diff(y))[same]), q01 = q[1], q025 = q[2],
                             stringsAsFactors = FALSE)
  }
  out <- list(seg = do.call(rbind, segs), bins = data.frame(chrom = chrom, start = start, end = end),
              stats = do.call(rbind, stats))
  if (!bin_noise) { out$seg$seg.sd <- NULL; out["stats"] <- list(NULL) }
  if (!is.null(prefix)) {
    wt <- function(d, f) utils::write.table(format(d, scientific = FALSE, trim = TRUE, digits = 15), f, sep = "\t",
                                            quote = FALSE, row.names = FALSE)
    wt(out$seg, paste0(prefix, ".seg")); wt(out$bins, paste0(prefix, ".bins.tsv"))
    if (bin_noise) wt(out$stats, paste0(prefix, ".sample_stats.tsv"))
  }
  out
}

#' Parameterise one tumour
#'
#' Steps 1 to 3 of the specification for a single sample: mode centring and the eight parameters
#' E and beta of the four detectors.
#'
#' @param mean,sd,nmark segment means (log2), standard deviations and marker counts.
#' @param dn,q01,q025 probe-level statistics; `NA` triggers the documented seg-only fallbacks.
#' @param params list from [rae_params()].
#' @param stream sub-stream of the random generator (the pipeline uses the 0-based sample index).
#' @return list with `E`, `beta`, `shift`, half-maximum points, the chosen peaks, flags and the density.
#' @export
rae_parameterize <- function(mean, sd = rep(NA_real_, length(mean)), nmark, dn = NA_real_, q01 = NA_real_,
                             q025 = NA_real_, params = rae_params(), stream = 0) {
  rae_parameterize_cpp(as.numeric(mean), as.numeric(sd), as.integer(nmark), dn, q01, q025, params, stream)
}

#' Unified breakpoint profile on one chromosome
#'
#' @param starts list with, per sample, the 1-based marker index at which each segment starts.
#' @param nmark list with the marker count of each of those segments.
#' @param n_markers number of markers on the chromosome.
#' @param params list from [rae_params()].
#' @return data frame of regions (`start`, `end`, 1-based inclusive marker indices).
#' @export
rae_ubp <- function(starts, nmark, n_markers, params = rae_params()) {
  rae_ubp_cpp(lapply(starts, as.integer), lapply(nmark, as.integer), as.integer(n_markers), params)
}

#' Regions of interest from region scores
#'
#' @param score,error,q summary score, its analytical error and the q-value of each region, in genomic order.
#' @param chrom integer chromosome index of each region.
#' @param params list from [rae_params()].
#' @return data frame with one row per stage-1 event and per stage-2 peak (1-based region indices).
#' @export
rae_roi <- function(score, error, q, chrom = rep(1L, length(score)), params = rae_params()) {
  rae_roi_cpp(as.numeric(score), as.numeric(error), as.numeric(q), as.integer(chrom), params)
}

#' @export
print.rae_result <- function(x, ...) {
  ev <- x$roi[x$roi$type == "region", , drop = FALSE]
  cat("RAE result:", nrow(x$samples), "samples,", nrow(x$regions), "regions in the unified breakpoint profile\n")
  cat("  significant events (q <=", x$params$fdr, "):", sum(ev$direction == "gain"), "gain,",
      sum(ev$direction == "loss"), "loss;", sum(x$roi$type == "peak"), "peaks\n")
  cat("  report:", grep("report.html$", x$files, value = TRUE), "\n")
  invisible(x)
}

#' Plot genome-wide significance
#'
#' Base-graphics counterpart of Figure 5 of the paper: -log10 q for gain (up) and loss (down).
#' The full figure set is in the HTML report listed in `x$files`.
#'
#' @param x a `rae_result`.
#' @param ... ignored.
#' @export
plot.rae_result <- function(x, ...) {
  r <- x$regions
  len <- tapply(r$end, r$chrom, max)
  off <- c(0, cumsum(as.numeric(len)))[match(r$chrom, names(len))]
  up <- -log10(r$A_q); dn <- log10(r$D_q)
  lim <- max(c(up, -dn, 2.5))
  graphics::plot(NA, xlim = c(0, sum(as.numeric(len))), ylim = c(-lim, lim), xaxt = "n", xlab = "chromosome",
                 ylab = "-log10 q   (loss below, gain above)", bty = "n")
  graphics::rect(off + r$start, 0, off + r$end, up, col = "#eb6834", border = NA)
  graphics::rect(off + r$start, 0, off + r$end, dn, col = "#2a78d6", border = NA)
  graphics::abline(h = c(-1, 1) * -log10(x$params$fdr), lty = 2, col = "#52514e")
  mids <- c(0, cumsum(as.numeric(len)))[seq_along(len)] + as.numeric(len) / 2
  graphics::axis(1, at = mids, labels = names(len), tick = FALSE, cex.axis = 0.7)
  graphics::abline(v = cumsum(as.numeric(len)), col = "#e1e0d9")
  invisible(x)
}
