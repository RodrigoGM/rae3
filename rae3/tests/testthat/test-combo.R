# All four opt-in analyses together (G-14a/b ploidy prior, G-14d focal_only, G-14e arm_test, G-14f
# peak_independence), on the liposarcoma cohort of helper-independence.R (hg38, so the arms are real arms).

combo_args <- list(arm_n_perm = 300L, arm_min_markers = 5L)
combo_run <- function(sf, out, ...)
  rae(read.delim(sf), build = "hg38", params = do.call(fast, c(combo_args, list(...))), out = out)
main6 <- c("regions", "roi", "samples", "segments", "calls", "null")

test_that("all four features: CLI = R; arm_test is the same with focal_only and leaves the main tables alone", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  skip_if_not_installed("jsonlite")
  dir <- tempfile("combo_"); dir.create(dir)
  sf <- file.path(dir, "c.seg")
  write.table(format(lipo_toy(), scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t", quote = FALSE,
              row.names = FALSE)
  rd <- function(p, t) readLines(file.path(dir, paste0(p, ".", t, ".tsv")))
  # (a) all four on: the seven tables are byte-identical, run.json equal apart from the input path
  all4 <- combo_run(sf, file.path(dir, "r"), ploidy_prior = TRUE, focal_only = TRUE, arm_test = TRUE,
                    peak_independence = TRUE)
  st <- system2(cli, c(sf, "--build", "hg38", "--n-perm", "400", "--recleave-every", "100", "--ploidy-prior",
                       "--focal-only", "--arm-test", "--arm-perm", "300", "--arm-min-markers", "5",
                       "--peak-independence", "--out", file.path(dir, "c")), stdout = FALSE, stderr = FALSE)
  expect_equal(st, 0L)
  for (t in c(main6, "arms")) expect_identical(rd("c", t), rd("r", t), info = t)
  js <- function(p) {
    j <- jsonlite::fromJSON(file.path(dir, paste0(p, ".run.json")), simplifyVector = FALSE)
    j$input <- NULL
    j
  }
  expect_identical(js("c"), js("r"))
  expect_true(any(!is.na(all4$roi$independent)))   # the independence test had a secondary peak to judge
  expect_gt(nrow(all4$arms), 0)
  # (b) the arm table does not depend on focal_only
  combo_run(sf, file.path(dir, "a"), arm_test = TRUE)
  combo_run(sf, file.path(dir, "af"), arm_test = TRUE, focal_only = TRUE)
  expect_identical(rd("a", "arms"), rd("af", "arms"))
  # (c) arm_test alone changes none of the six main tables (ruling R12)
  combo_run(sf, file.path(dir, "d"))
  for (t in main6) expect_identical(rd("a", t), rd("d", t), info = t)
})
