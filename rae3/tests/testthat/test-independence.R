# Independence of a secondary peak (gap register G-14f): inside one event, each later stage-2 peak is re-scored
# with the samples' runs through the accepted peaks' apexes masked, and called independent when the residual still
# reaches the significance cut. Opt-in; off changes nothing (test-defaults.R).

# ind_toy(), shoulder_toy() and lipo_toy() are in helper-independence.R.
ind_run <- function(d, preset = "sequencing", ...)
  rae(d, build = "hg38", params = fast(preset = preset, peak_independence = TRUE, ...), out = tempfile("ind_"))
gain_peaks <- function(res) { r <- res$roi; r[r$direction == "gain" & r$type == "peak", ] }

test_that("shoulder peak is not independent", {
  res <- ind_run(shoulder_toy())
  pk <- gain_peaks(res)
  pk <- pk[pk$chrom == 3, ]
  expect_equal(nrow(pk), 2)
  expect_equal(length(unique(pk$event)), 1)
  first <- which.max(pk$score)
  expect_true(is.na(pk$independent[first]))
  expect_true(is.na(pk$resid_score[first]) && is.na(pk$resid_p[first]))
  expect_identical(pk$independent[-first], 0L)
  # the same samples' runs carry both maxima: masking them leaves only A1 below the run threshold (0.28 x 12/40)
  cut <- min(res$regions$A_score[res$regions$A_q <= res$params$peak_q])
  expect_lt(pk$resid_score[-first], cut)
  expect_lt(pk$resid_score[-first], pk$score[-first] / 3)
  expect_gt(pk$resid_p[-first], 0.01)
})

test_that("two drivers carried by different samples are independent", {
  res <- ind_run(lipo_toy(), focal_only = TRUE)
  pk <- gain_peaks(res)
  pk <- pk[pk$chrom == 3, ]
  expect_equal(nrow(pk), 2)
  expect_equal(length(unique(pk$event)), 1)
  first <- which.max(pk$score)
  expect_true(is.na(pk$independent[first]))
  expect_identical(pk$independent[-first], 1L)
  # no run of the second driver's samples reaches the first apex: nothing is masked, so the residual is the score
  expect_equal(pk$resid_score[-first], pk$score[-first], tolerance = 1e-5)
  expect_lte(pk$resid_p[-first], 0.01)
  # without focal_only the shared arm gain is one run through both apexes, so the second peak is not independent
  d <- gain_peaks(ind_run(lipo_toy()))
  d <- d[d$chrom == 3, ]
  expect_equal(nrow(d), 2)
  expect_identical(d$independent[-which.max(d$score)], 0L)
})

test_that("primary and single peaks are NA", {
  # one amplicon at 100-110 Mb with a single maximum at 103-105 Mb: an event of several regions, one stage-2 peak.
  # Every preset, including cgh, whose stage-2 peaks are the coded detector's (roi_code).
  one <- ind_toy(ev = list(list(s = 1:12, a = 100e6, b = 110e6, v = 0.6), list(s = 1:12, a = 103e6, b = 105e6, v = 0.8)))
  for (pr in c("sequencing", "cgh", "taylor2008")) {
    r <- ind_run(one, preset = pr)$roi
    expect_true(all(c("resid_score", "resid_p", "independent") %in% names(r)), info = pr)
    expect_type(r$independent, "integer")         # typed although every value is NA
    expect_type(r$resid_score, "double")
    expect_type(r$resid_p, "double")
    pk <- r[r$type == "peak" & r$chrom == 3, ]
    expect_equal(nrow(pk), 1, info = pr)
    expect_gt(r$n_regions[r$type == "region" & r$event == pk$event & r$direction == pk$direction], 1)
    expect_true(all(is.na(r$independent) & is.na(r$resid_score) & is.na(r$resid_p)), info = pr)
  }
  # cgh, coded peaks on the shoulder cohort: the columns are written and the primary of every event is NA
  g <- ind_run(shoulder_toy(), preset = "cgh")$roi
  expect_true(all(c("resid_score", "resid_p", "independent") %in% names(g)))
  pk <- g[g$type == "peak", ]
  expect_gt(nrow(pk), 0)
  for (k in split(pk, paste(pk$direction, pk$event))) expect_true(is.na(k$independent[which.max(k$score)]))
  expect_true(all(is.na(g$independent[g$type == "region"])))
  # the primary is the highest peak, not the first in genomic order: here the second peak is higher
  two <- ind_toy(ev = list(list(s = 1:12, a = 95e6, b = 115e6, v = 0.4), list(s = 1:12, a = 98e6, b = 100e6, v = 1.2),
                           list(s = 13:20, a = 110e6, b = 112e6, v = 1.6)))
  pk <- gain_peaks(ind_run(two))
  pk <- pk[pk$chrom == 3, ]
  expect_equal(nrow(pk), 2)
  expect_gt(pk$score[2], pk$score[1])
  expect_true(is.na(pk$independent[2]) && is.na(pk$resid_score[2]))
  # samples 1-12: the A0 run of their 95-115 Mb gain reaches the primary apex and is masked; their A1 run (the
  # amplification at 98-100 Mb) does not, so the peak keeps its A1 part and stays above the cut
  expect_identical(pk$independent[1], 1L)
  expect_lt(pk$resid_score[1], pk$score[1])
})

test_that("columns only when on", {
  hdr <- function(res) readLines(grep("roi\\.tsv$", res$files, value = TRUE), n = 1)
  off <- rae(toy_cohort(n = 20), build = "custom", params = fast(seed = 1L), out = tempfile("ind_"))
  expect_identical(hdr(off), readLines(testthat::test_path("fixtures", "defaults", "seq_custom", "roi.tsv"), n = 1))
  expect_false(any(c("resid_score", "resid_p", "independent") %in% names(off$roi)))
  on <- rae(toy_cohort(n = 20), build = "custom", params = fast(seed = 1L, peak_independence = TRUE), out = tempfile("ind_"))
  expect_identical(hdr(on), paste0(hdr(off), "\tresid_score\tresid_p\tindependent"))
  # everything else in roi.tsv is unchanged
  strip <- function(res) sub("(\t[^\t]*){3}$", "", readLines(grep("roi\\.tsv$", res$files, value = TRUE))[-1])
  expect_identical(strip(on), readLines(grep("roi\\.tsv$", off$files, value = TRUE))[-1])
  expect_identical(on$params$peak_independence, TRUE)
  expect_identical(off$params$peak_independence, FALSE)
  expect_identical(off$params$peak_run_min, 0.5)
})

test_that("peak_run_min is validated", {
  expect_error(ind_run(toy_cohort(n = 10), peak_run_min = 0), "peak_run_min")
  expect_error(ind_run(toy_cohort(n = 10), peak_run_min = 1.5), "peak_run_min")
})

test_that("CLI and R package agree with --peak-independence --peak-run-min 0.4", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  dir <- tempfile("iapp_"); dir.create(dir)
  sf <- file.path(dir, "c.seg")
  write.table(format(lipo_toy(), scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t",
              quote = FALSE, row.names = FALSE)
  rae(read.delim(sf), build = "hg38", params = fast(focal_only = TRUE, peak_independence = TRUE, peak_run_min = 0.4),
      out = file.path(dir, "r.rae"))
  st <- system2(cli, c(sf, "--build", "hg38", "--n-perm", "400", "--recleave-every", "100", "--focal-only",
                       "--peak-independence", "--peak-run-min", "0.4", "--out", file.path(dir, "c.rae")),
                stdout = FALSE, stderr = FALSE)
  expect_equal(st, 0)
  for (t in c("regions", "roi", "samples", "segments", "calls", "null"))
    expect_identical(readLines(file.path(dir, paste0("c.rae.", t, ".tsv"))),
                     readLines(file.path(dir, paste0("r.rae.", t, ".tsv"))), info = t)
  cr <- rae_read(file.path(dir, "c.rae"))
  expect_identical(cr$params$peak_run_min, 0.4)
  expect_true(any(cr$roi$independent %in% 1L))
  expect_true(any(grepl("--peak-independence", system2(cli, "--help", stdout = TRUE))))
})

test_that("the masked recompute with nothing masked gives the observed region scores", {
  # sequencing (marker weights), cgh (D-14 squared weights), and the focal scores
  for (p in list(fast(), fast(preset = "cgh"), fast(focal_only = TRUE))) {
    m <- rae3:::rae_masked_means_cpp(rae3:::.seg_df(lipo_toy()), NULL, NULL, "hg38", p)
    expect_gt(sum(!is.na(m$observed)), 0)
    expect_identical(m$recomputed, m$observed)
  }
})
