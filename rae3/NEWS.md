# rae3 (development version)

## New features

* Renamed `rae2` → `rae3` (version 3.0.0.9000): RAE beta 0.91 (2008) and Rae V2 (N. Socci, 2010–2015)
  came first.
* Presets: `sequencing` (the default; see below) and two for array CGH. `cgh` is the released RAE code
  (beta 0.91), which reproduces the published results; `taylor2008` is a clean-room reading of the 2008
  paper's text. The 15 places where that reading and the code differ are details the text leaves open,
  documented in `03_original_code_diff.md`.
* Noise can be measured from QDNAseq bins: `bin_noise = TRUE` in `rae()` and `as_seg()`. The default is
  `FALSE`, so R and the command line agree on plain `.seg` input.
* Regions are annotated for overlap with common population CNVs: DGV Gold Standard (hg19), gnomAD-SV v4.1
  (hg38), and HPRC long-read deletions (T2T-CHM13).
* `filter_germline_overlaps`: when `FALSE` (the default), regions that overlap germline variation are
  drawn grey and marked "Г". When `TRUE`, they move to `roi_germline`.
* Plotting functions: `rae_plot_genome()`, `rae_plot_locus()` and `rae_heatmap()`, a noise plot, and the
  diagnostic figures of Taylor et al. 2008. Save them with `rae_save()` and style them with
  `theme_rodrigo()`.
* Cancer-gene labels from OncoKB, fetched locally: `rae_fetch_cancer_genes()`.
* Preset `sequencing`: per-tumour detectors placed on the copy-number ladder instead of from the noise
  model, which fails on sequencing data (70-80 % of BriTROC-1 tumours fall back to fixed slopes). The one-copy
  step comes from supplied purity/ploidy (`purity_ploidy` in `rae()`, `--purity-ploidy` on the command line)
  or is estimated from the periodicity of the segment means. On BriTROC-1 the amplification agreement with the
  published calls (kappa) is 0.65 estimated and 0.95 with purity/ploidy, against 0.52 for `cgh`.
* `rae_plot_frequency()`: percentage of tumours with each alteration, stacked or one state at a time.
* `rae_read()` loads the command-line tool's output into a `rae()` result, so every plot works on
  command-line runs. The output set gains an eighth file, `.null.tsv`, with the null-model internals.
* Two fractions of the genome altered (FGA) per tumour in `samples`: `rae_fga`, by the RAE calls (the
  Taylor 2010 / Hieronymus 2014 CNA burden), and `fga`, by a fixed threshold, |log2 ratio| > `fga_threshold`
  (0.2, the cBioPortal / MSK-IMPACT FGA; `fga_threshold_loss` sets the loss side separately, e.g. for
  single-cell data). Also `n_altered`, the number of segments counted in `rae_fga`.
  See *Output* in `?rae`.

### Opt-in analyses

All four are off by default and change no default output (`test-defaults.R` compares byte for byte). Each is a
switch in `rae_params()` and a flag on the command line; see `?rae_params`, *Output* in `?rae`, and gap
register G-14a to G-14f.

* `focal_only` (`--focal-only`): each segment is scored on its deviation from its arm's baseline, with
  segments cut at the centromere, so arm-level change no longer reaches the scores, the null model or the
  regions. Runs of deviation longer than `focal_max_bp` (20 Mb) score 0, and `samples` gains
  `focal_excluded_frac`. `rae_fga`, `fga` and `n_altered` stay absolute. The arm baseline is set by
  `arm_baseline`: `"median"` (default, the marker-weighted median of the arm's segments) or `"twolevel"`
  (a one- or two-level fit per chromosome); `arm_min_markers` (50) sets which arms are eligible. Median was
  chosen by a pre-registered comparison, but both truths used (BriTROC-1 absolute copy number, a simulation)
  are arm-majority states, which the median estimates directly; a fairer test of `twolevel` is planned.
* `arm_test` (`--arm-test`): a cohort test of each chromosome arm by within-tumour permutation
  (`arm_n_perm`, 10000), with permutation p, Benjamini-Hochberg q and Westfall-Young max-T p for gain and
  loss. It writes a ninth file, `<prefix>.arms.tsv`, only when on; `rae()` returns it as `arms` and
  `rae_plot_arms()` draws it. `rae_read()` reads the file only when `run.json` says the test was on.
  Where the build has no centromere, each chromosome is one arm and its label is the bare chromosome number.
* `peak_independence` (`--peak-independence`, with `peak_run_min`, 0.5): for events with more than one peak,
  each later peak is re-scored with the alterations through the higher peaks removed. `roi` gains
  `resid_score`, `resid_p` (a heuristic comparison with the observed cut, not a calibrated p-value) and
  `independent`. Stretches below `peak_run_min` are never removed, so a weak shoulder can still be called
  independent.
* `ploidy_prior` (`--ploidy-prior`, with `ploidy_prior_wgd`, `_mu1`, `_mu2`, `_sd`, `_ambig`): a soft prior on
  mean ploidy (two normal components, diploid-like and doubled genome) chooses the ladder's base copy number
  instead of the lowest-level rule. `samples` gains `ploidy_hard`, `ploidy_alt`, `ploidy_prior_ratio` (NA
  for tumours with a single admissible candidate) and the flag `ploidy_ambiguous`. On BriTROC-1 the
  pre-registered benchmark found no improvement: at the defaults it changed no tumour, and 98 of 265
  tumours are flagged ambiguous (`docs/spec/08_britroc.md`, section 8).

## Changes

* **The default preset is now `sequencing`**, so results change for calls that pass no preset. For array
  data pass `preset = "cgh"` (R) or `--preset cgh` (command line).
* Presets renamed during development: `paper` → `taylor2008`, `original` → `cgh`. The old names are
  rejected; `rae_read()` still reads runs saved under them.
* **`cna_burden` is renamed `rae_fga`** (a fraction of the genome altered, by the RAE calls) **and now
  divides by the tumour's segment-covered span** (first to last segment on each autosome), the published
  definition, instead of the whole assembly. Values rise by about 3 % on arrays (×1.03) under every
  preset; the ratio to Hieronymus 2014 and Taylor 2010 goes from 0.967 to 0.998, and 145 instead of 141 of
  147 Taylor 2010 tumours fall on the published side of the 5.41 % cut-off. `fga_assembly = TRUE`
  (`--fga-assembly`) restores the old denominator; `rae_read()` maps old runs. See *Output* in `?rae`.
* `--preset` is applied before the other command-line options wherever it appears, so an override such as
  `--n-perm` is never reset by a later preset.
* BSD 3-Clause licence.
* hg17 is no longer a user-facing build; it is kept for development only.
* Copy number is reported as log2 ratio throughout.
* `samples.tsv` gains `copy_step_source`, `copy_step`, `purity`, `ploidy`; a new flag `step_fallback`.
* `run.json` is valid JSON (missing values are `null`) and records parameters at full precision.

## Validation

* The published per-tumour CNA burden is reproduced: Pearson 0.9999 for Hieronymus 2014 and 0.997 for
  Taylor 2010 (`cgh` preset).
* Taylor 2010 focal regions (Table S3): 58 of 66 recovered. Arm-level events (Table S4): 7 of 7.

# rae2 2.0.0

* First implementation of RAE from the publication: a C++ core shared by the R package and the `rae`
  command-line tool. It takes QDNAseq objects and `.seg` input.
