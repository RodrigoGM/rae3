# Cases for the "defaults unchanged" regression test (test-defaults.R) and for the script that
# produced its fixtures (tools/make_default_fixtures.R). Needs helper-toy.R (toy_cohort, fast).

default_cases <- function() list(
  seq_custom   = list(preset = "sequencing", build = "custom"),
  cgh_custom   = list(preset = "cgh",        build = "custom"),
  t2008_custom = list(preset = "taylor2008", build = "custom"),
  seq_hg38     = list(preset = "sequencing", build = "hg38"))

# Runs one case on the toy cohort (seed 1) and returns the directory holding toy.*.tsv files.
run_case <- function(case, ...) {
  dir <- tempfile("rae_case_"); dir.create(dir)
  rae(toy_cohort(n = 20), build = case$build, params = fast(preset = case$preset, seed = 1L, ...),
      out = file.path(dir, "toy"))
  dir
}
