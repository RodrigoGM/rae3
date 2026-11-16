# Arm geometry, centromere split and the two arm-baseline estimators (G-14c). Internal hooks, called as rae3:::.

# 1 Mb bins over all of chr1 (hg38: centromere 121.7-125.1 Mb) and the first 20 Mb of chr2
arm_bins <- function() {
  data.frame(chrom = c(rep("1", 248), rep("2", 20)),
             start = c(0:247 * 1e6 + 1, 0:19 * 1e6 + 1), end = c(1:248 * 1e6, 1:20 * 1e6))
}
hook <- function(values, weights, chrom, arm, method) {
  rae3:::rae_arm_baseline_cpp(values, weights, chrom, arm, method)
}

test_that("arm layout splits at the centromere and marks eligibility", {
  a <- rae3:::rae_arm_layout_cpp("hg38", 124L, arm_bins())
  a1 <- a[a$chrom == 1, ]
  expect_identical(a1$arm, c("p", "q"))
  # marker midpoints 0.5, 1.5, ... Mb: < 121.7 is p; 122.5 (inside) is nearer the start; 123.5 nearer the end
  expect_identical(a1$m_first, c(1L, 124L))
  expect_identical(a1$m_last, c(123L, 248L))
  expect_identical(a1$n_markers, c(123L, 125L))
  expect_identical(a1$eligible, c(FALSE, TRUE))          # needs >= 124 markers
  expect_equal(a1$bp_start, c(1, 123e6 + 1))
  expect_equal(a1$bp_end, c(123e6, 248e6))
  expect_true(all(rae3:::rae_arm_layout_cpp("hg38", 50L, arm_bins())$eligible[1:2]))
  # chr2 bins end at 20 Mb, all left of the centromere: a single p arm; markers continue the global index
  a2 <- a[a$chrom == 2, ]
  expect_identical(a2$arm, "p")
  expect_identical(a2$m_first, 249L)
  # unknown centromere: one whole-chromosome arm
  w <- rae3:::rae_arm_layout_cpp("custom", 50L, arm_bins())
  expect_identical(w$arm, c("w", "w"))
  expect_identical(w$n_markers, c(248L, 20L))
})

test_that("segments crossing the centromere are split with nmark conserved", {
  seg <- data.frame(ID = "S1", chrom = "1",
                    loc.start = c(1, 100e6 + 1, 130e6 + 1), loc.end = c(100e6, 130e6, 248e6),
                    num.mark = c(100L, 30L, 10L), seg.mean = c(0, 0.5, -0.2))
  o <- rae3:::rae_split_segments_cpp(seg, arm_bins()[1:248, ], "hg38", 50L)
  expect_equal(nrow(o), 4)
  expect_equal(sum(o$num.mark), sum(seg$num.mark))
  # the middle segment spans markers 101-130: 101-123 are p (23), 124-130 q (7); nmark 30 splits 23 : 7
  mid <- o[o$seg.mean == 0.5, ]
  expect_equal(mid$num.mark, c(23, 7))
  expect_equal(mid$loc.start, c(100e6 + 1, 123e6 + 1))
  expect_equal(mid$loc.end, c(123e6, 130e6))
  # a segment with few markers: each side still gets at least one, and the total is kept
  seg2 <- data.frame(ID = "S1", chrom = "1", loc.start = 120e6 + 1, loc.end = 128e6, num.mark = 3L, seg.mean = 1)
  o2 <- rae3:::rae_split_segments_cpp(seg2, arm_bins()[1:248, ], "hg38", 50L)
  expect_equal(nrow(o2), 2)
  expect_true(all(o2$num.mark >= 1)); expect_equal(sum(o2$num.mark), 3)
  # a segment with fewer than 2 markers cannot become two real pieces: it is left whole (judged on the arm of its start marker)
  seg3 <- data.frame(ID = "S1", chrom = "1", loc.start = 120e6 + 1, loc.end = 128e6, num.mark = 1L, seg.mean = 1)
  o3 <- rae3:::rae_split_segments_cpp(seg3, arm_bins()[1:248, ], "hg38", 50L)
  expect_equal(nrow(o3), 1)
  expect_equal(o3$num.mark, 1); expect_equal(o3$loc.start, 120e6 + 1); expect_equal(o3$loc.end, 128e6)
})

test_that("the centromere cut keeps the segment's bp length when the grid has a gap there", {
  # no markers over 121-126 Mb (around the hg38 chr1 centromere, 121.7-125.1 Mb)
  keep <- c(0:120, 126:247)
  gb <- data.frame(chrom = "1", start = keep * 1e6 + 1, end = (keep + 1) * 1e6)
  seg <- data.frame(ID = "S1", chrom = "1", loc.start = c(1, 100e6 + 1, 130e6 + 1), loc.end = c(100e6, 130e6, 248e6),
                    num.mark = c(100L, 25L, 118L), seg.mean = c(0, 0.5, -0.2))
  o <- rae3:::rae_split_segments_cpp(seg, gb, "hg38", 50L)
  mid <- o[o$seg.mean == 0.5, ]
  expect_equal(nrow(mid), 2)
  expect_equal(sum(mid$loc.end - mid$loc.start + 1), 30e6)
  expect_equal(mid$loc.start, c(100e6 + 1, 126e6 + 1))       # cut_bp = first q marker's start - 1
  expect_equal(mid$loc.end, c(126e6, 130e6))
  expect_equal(sum(o$loc.end - o$loc.start + 1), sum(seg$loc.end - seg$loc.start + 1))
})

test_that("median baseline equals the weighted median", {
  # lower weighted median: smallest value whose cumulative weight reaches half the total
  o <- hook(c(0, 1, 2), c(1, 1, 5), c(1L, 1L, 1L), c("p", "p", "p"), 0L)
  expect_equal(o$level, 2)
  o <- hook(c(0, 1), c(1, 1), c(1L, 1L), c("p", "p"), 0L)
  expect_equal(o$level, 0)
  o <- hook(c(3, 0, 1, 2, 9), c(1, 1, 1, 1, 1), rep(1L, 5), rep("q", 5), 0L)
  expect_equal(o$level, 2)
  # one level per arm present in the data
  o <- hook(c(0, 0, 5), c(1, 1, 1), c(1L, 1L, 1L), c("p", "p", "q"), 0L)
  expect_equal(o$arm, c("p", "q")); expect_equal(o$level, c(0, 5))
})

test_that("two-level fit finds a mid-chromosome step and keeps one level on a flat chromosome", {
  v <- c(rep(0, 6), rep(1, 14))
  arm <- c(rep("p", 10), rep("q", 10))
  o <- hook(v, rep(10, 20), rep(1L, 20), arm, 1L)
  expect_equal(o$level, c(0, 1))                      # p holds 6 segments at 0 and 4 at 1: the 0 level covers more
  expect_equal(o$coverage, c(0.6, 1))
  # the arm median is the same here; the two methods differ only when the step is not at the arm's majority
  expect_equal(hook(v, rep(10, 20), rep(1L, 20), arm, 0L)$level, c(0, 1))
  # a flat chromosome: both arms share one level
  flat <- hook(rep(0.3, 20), rep(10, 20), rep(1L, 20), arm, 1L)
  expect_equal(flat$level, c(0.3, 0.3)); expect_equal(flat$coverage, c(1, 1))
  set.seed(5); n <- rnorm(20, 0, 0.01)
  fl <- hook(n, rep(10, 20), rep(1L, 20), arm, 1L)
  expect_equal(fl$level[1], fl$level[2])
  # chromosomes are fitted separately; a one-segment chromosome is a single level
  two <- hook(c(0, 0, 2), c(5, 5, 5), c(1L, 1L, 2L), c("p", "q", "p"), 1L)
  expect_equal(two$level, c(0, 0, 2))
  # a step between two segments is found even when the arm medians would tie on weights
  st <- hook(c(0, 0, 0, 3, 3, 3), c(100, 1, 1, 1, 1, 100), rep(1L, 6), c(rep("p", 3), rep("q", 3)), 1L)
  expect_equal(st$level, c(0, 3))
})

test_that("arm_baseline is a character argument stored as the core's code and round-trips through the manifest", {
  expect_identical(rae_params()$arm_baseline, 0L)   # median by default
  expect_identical(rae_params(arm_baseline = "twolevel")$arm_baseline, 1L)
  expect_identical(rae_params(arm_baseline = "median")$arm_baseline, 0L)
  expect_error(rae_params(arm_baseline = "mean"))
  res <- rae(toy_cohort(), build = "hg19", params = fast(arm_baseline = "twolevel", arm_min_markers = 20L),
             out = tempfile("arm_"))
  rd <- rae_read(sub("\\.regions\\.tsv$", "", grep("regions", res$files, value = TRUE)))
  expect_equal(rd$params$arm_baseline, 1L)
  expect_equal(rd$params$arm_min_markers, 20L)
  expect_equal(rd$params, res$params)
})

test_that("baseline coverage is 1 on a flat arm", {
  o <- hook(rep(2, 8), rep(3, 8), rep(1L, 8), rep("p", 8), 1L)
  expect_equal(o$coverage, 1)
  o <- hook(rep(2, 8), rep(3, 8), rep(1L, 8), rep("p", 8), 0L)
  expect_equal(o$coverage, 1)
})

test_that("arm baselines use copy units with a ladder and the detectors without one", {
  b <- arm_bins()[1:248, ]
  prm <- list(arm_min_markers = 1L, arm_baseline = 0L)
  # the first segment (markers 1-200) is cut at marker 123; its q piece (77 markers) meets a 3-copy segment (48)
  seg <- data.frame(ID = "S1", chrom = "1", loc.start = c(1, 200e6 + 1), loc.end = c(200e6, 248e6),
                    num.mark = c(200L, 48L), seg.mean = c(0, log2(1.5)))
  # ladder: ratio = m * (r0 + c * step); m = 1, r0 = 0, step = 0.5: ratio 1 is 2 copies, 1.5 is 3 copies
  o <- rae3:::rae_arm_baselines_cpp(seg, b, "hg38", prm, step = 0.5, lad_m = 1, lad_r0 = 0, E = numeric(), B = numeric())
  expect_equal(o$arm, c("p", "q"))
  expect_equal(o$level, c(2, 2))
  expect_equal(o$coverage, c(1, 200 * 77 / 200 / (200 * 77 / 200 + 48)))
  # log2 scale with detectors at +-0.3: the 3-copy piece (log2 1.5 = 0.58) is above E_A0, so not at baseline
  E <- c(0.3, 0.6, -0.3, -0.6); B <- c(20, 20, -20, -20)
  o <- rae3:::rae_arm_baselines_cpp(seg, b, "hg38", prm, step = NA, lad_m = NA, lad_r0 = NA, E = E, B = B)
  expect_equal(o$level, c(0, 0))
  expect_equal(o$coverage, c(1, 77 / 125))
})

# ---- arm-level cohort test (G-14e) ------------------------------------------------------------------------
# hg38 chromosomes 1-6 at full length, 1 Mb markers, one segment per arm cut at a whole Mb inside the
# centromere. `alt[[j]]` is a named vector arm -> log2 value for sample j (other arms 0); `drop[[j]]` lists arms
# with no segment in sample j.
armtest_genome <- function(chroms = 1:6) {
  cb <- rae_cytobands("hg38")
  data.frame(chrom = chroms,
             cen = vapply(chroms, function(c) round(cb$end[cb$chrom == c & cb$stain == "acen"][1] / 1e6) * 1e6, 0),
             len = vapply(chroms, function(c) floor(max(cb$end[cb$chrom == c]) / 1e6) * 1e6, 0))
}
armtest_toy <- function(alt, drop = list(), chroms = 1:6, seed = 1, bs = 1e6) {
  set.seed(seed)
  g <- armtest_genome(chroms)
  rows <- list()
  for (j in seq_along(alt)) for (i in seq_len(nrow(g))) for (a in c("p", "q")) {
    lab <- paste0(g$chrom[i], a)
    if (j <= length(drop) && lab %in% drop[[j]]) next
    lo <- if (a == "p") 0 else g$cen[i]
    hi <- if (a == "p") g$cen[i] else g$len[i]
    v <- if (lab %in% names(alt[[j]])) alt[[j]][[lab]] else 0
    rows[[length(rows) + 1]] <- data.frame(ID = sprintf("A%02d", j), chrom = g$chrom[i], loc.start = lo + 1, loc.end = hi,
                                           num.mark = round((hi - lo) / bs), seg.mean = v + rnorm(1, 0, 0.01), seg.sd = 0.12)
  }
  do.call(rbind, rows)
}
armtest_labels <- function(chroms = 1:6) paste0(rep(chroms, each = 2), c("p", "q"))
# one random arm per sample (excluding `avoid`), gain +0.5 or loss -0.6
armtest_passengers <- function(n, avoid = character(0), chroms = 1:6)
  lapply(seq_len(n), function(j) stats::setNames(sample(c(0.5, -0.6), 1), sample(setdiff(armtest_labels(chroms), avoid), 1)))
armtest_run <- function(d, out = tempfile("armt_"), ...) {
  a <- utils::modifyList(list(arm_test = TRUE, arm_n_perm = 500L, arm_min_markers = 20L), list(...))
  rae(d, build = "hg38", params = do.call(fast, a), out = out)
}

test_that("planted arm loss is significant, unaltered arms are not", {
  set.seed(21)
  alt <- armtest_passengers(30, avoid = "2p")
  for (j in 1:15) alt[[j]] <- c(alt[[j]], `2p` = -0.6)
  a <- armtest_run(armtest_toy(alt, seed = 21))$arms
  expect_equal(nrow(a), 12)
  expect_identical(a$arm_label, armtest_labels())
  expect_lte(a$q_loss[a$arm_label == "2p"], 0.01)
  expect_lte(a$fwer_loss[a$arm_label == "2p"], 0.01)
  expect_gt(median(a$p_loss[a$arm_label != "2p"]), 0.2)
  expect_equal(a$frac_loss[a$arm_label == "2p"], 0.5)
  expect_true(all(a$n_samples == 30))
})

test_that("arm p-values are calibrated under random arm placement", {
  # 20 cohorts of 16 tumours, each with one altered arm chosen at random: every arm is exchangeable, so the
  # pooled-null p is uniform. All arms of all cohorts are pooled (within a cohort the p-values are negatively
  # dependent, as ranks are, which makes the KS test conservative, not liberal).
  pg <- pl <- numeric(0)
  for (s in 1:20) {
    set.seed(100 + s)
    a <- armtest_run(armtest_toy(armtest_passengers(16), seed = 100 + s), seed = s)$arms
    pg <- c(pg, a$p_gain); pl <- c(pl, a$p_loss)
  }
  expect_gt(suppressWarnings(stats::ks.test(pg, "punif"))$p.value, 0.01)
  expect_gt(suppressWarnings(stats::ks.test(pl, "punif"))$p.value, 0.01)
})

test_that("fwer_gain >= p_gain and fwer is monotone in score", {
  set.seed(5)
  alt <- armtest_passengers(24)
  for (j in 1:8) alt[[j]] <- c(alt[[j]], `3q` = 0.5)
  a <- armtest_run(armtest_toy(alt, seed = 5))$arms
  expect_true(all(a$fwer_gain >= a$p_gain))
  expect_true(all(a$fwer_loss >= a$p_loss))
  expect_true(all(a$q_gain >= a$p_gain - 1e-12))
  for (d in c("gain", "loss")) {
    o <- order(a[[paste0("score_", d)]])
    expect_true(all(diff(a[[paste0("fwer_", d)]][o]) <= 0), info = d)
    expect_true(all(diff(a[[paste0("p_", d)]][o]) <= 0), info = d)
  }
  expect_equal(a$score_gain, sqrt(a$A0^2 + a$A1^2), tolerance = 1e-5)
  expect_equal(a$score_loss, sqrt(a$D0^2 + a$D1^2), tolerance = 1e-5)
})

test_that("missing arm data are skipped, not zero", {
  set.seed(8)
  alt <- armtest_passengers(30, avoid = "1q")
  for (j in 1:10) alt[[j]] <- c(alt[[j]], `1q` = 0.5)
  drop <- c(rep(list(character(0)), 25), rep(list("1q"), 5))
  d <- armtest_toy(alt, drop = drop, seed = 8)
  a <- armtest_run(d)$arms
  q1 <- a[a$arm_label == "1q", ]
  expect_equal(q1$n_samples, 25)
  expect_true(all(a$n_samples[a$arm_label != "1q"] == 30))
  # the same 25 tumours alone (same order, so the same per-tumour streams) give the same 1q means
  b <- armtest_run(d[d$ID %in% sprintf("A%02d", 1:25), ])$arms
  b1 <- b[b$arm_label == "1q", ]
  for (k in c("A0", "A1", "D0", "D1", "frac_gain", "baseline_coverage")) expect_equal(q1[[k]], b1[[k]], info = k)
  expect_equal(q1$frac_gain, 10 / 25)
})

test_that("arms.tsv only when arm_test; rae_read returns x$arms; main analysis unchanged", {
  set.seed(9)
  d <- armtest_toy(armtest_passengers(16), seed = 9)
  dir <- tempfile("armio_"); dir.create(dir)
  on <- armtest_run(d, out = file.path(dir, "on"))
  off <- rae(d, build = "hg38", params = fast(arm_min_markers = 20L), out = file.path(dir, "off"))
  expect_true(file.exists(file.path(dir, "on.arms.tsv")))
  expect_false(file.exists(file.path(dir, "off.arms.tsv")))
  expect_null(off$arms)
  expect_s3_class(on$arms, "data.frame")
  expect_identical(names(on$arms), c("arm", "chrom", "arm_label", "start", "end", "length_bp", "n_markers", "n_samples",
    "A0", "A1", "D0", "D1", "score_gain", "score_loss", "p_gain", "q_gain", "fwer_gain", "p_loss", "q_loss",
    "fwer_loss", "frac_gain", "frac_loss", "baseline_coverage"))
  expect_true(file.path(dir, "on.arms.tsv") %in% on$files)
  rd <- rae_read(file.path(dir, "on"))
  expect_identical(rd$arms, on$arms)
  expect_true(rd$params$arm_test); expect_identical(rd$params$arm_n_perm, 500L)
  # the six data tables are byte-identical to a run without the arm test
  for (t in c("regions", "roi", "samples", "segments", "calls", "null"))
    expect_identical(readLines(file.path(dir, paste0("on.", t, ".tsv"))), readLines(file.path(dir, paste0("off.", t, ".tsv"))), info = t)
  skip_if_not_installed("jsonlite")
  js <- jsonlite::fromJSON(file.path(dir, "on.run.json"))
  expect_identical(js$summary$arms_significant_gain, sum(on$arms$q_gain <= 0.01))
  expect_identical(js$summary$arms_significant_loss, sum(on$arms$q_loss <= 0.01))
  expect_null(jsonlite::fromJSON(file.path(dir, "off.run.json"))$summary$arms_significant_gain)
  # with focal_only on as well the arm test reuses the same baselines: the same table
  fo <- armtest_run(d, focal_only = TRUE)
  expect_identical(fo$arms, on$arms)
})

test_that("a run without arm_test removes an earlier arms.tsv under the same prefix", {
  set.seed(11)
  d <- armtest_toy(armtest_passengers(10), seed = 11)
  out <- file.path(tempfile("armstale_"), "x")
  expect_false(is.null(armtest_run(d, out = out)$arms))
  expect_true(file.exists(paste0(out, ".arms.tsv")))
  off <- rae(d, build = "hg38", params = fast(arm_min_markers = 20L), out = out)
  expect_null(off$arms)
  expect_false(file.exists(paste0(out, ".arms.tsv")))
  expect_false(any(grepl("arms\\.tsv$", off$files)))
  # rae_read goes by run.json: a stray arms.tsv beside an arm_test-off run is not read
  writeLines("arm\tchrom", paste0(out, ".arms.tsv"))
  rd <- rae_read(out)
  expect_null(rd$arms)
  expect_false(any(grepl("arms\\.tsv$", rd$files)))
})

test_that("the CLI removes an earlier arms.tsv when run without --arm-test", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  set.seed(13)
  dir <- tempfile("armclistale_"); dir.create(dir)
  sf <- file.path(dir, "c.seg"); out <- file.path(dir, "c.rae")
  write.table(format(armtest_toy(armtest_passengers(10), seed = 13), scientific = FALSE, trim = TRUE, digits = 15), sf,
              sep = "\t", quote = FALSE, row.names = FALSE)
  base <- c(sf, "--build", "hg38", "--n-perm", "400", "--recleave-every", "100", "--arm-min-markers", "20", "--out", out)
  system2(cli, c(base, "--arm-test", "--arm-perm", "100"), stdout = FALSE, stderr = FALSE)
  expect_true(file.exists(paste0(out, ".arms.tsv")))
  system2(cli, base, stdout = FALSE, stderr = FALSE)
  expect_false(file.exists(paste0(out, ".arms.tsv")))
  expect_true(file.exists(paste0(out, ".regions.tsv")))
  expect_null(rae_read(out)$arms)
})

test_that("arm seed is exposed", {
  set.seed(10)
  d <- armtest_toy(armtest_passengers(16), seed = 10)
  dir <- tempfile("armseed_"); dir.create(dir)
  armtest_run(d, out = file.path(dir, "a")); armtest_run(d, out = file.path(dir, "b"))
  armtest_run(d, out = file.path(dir, "c"), seed = 2)
  rd <- function(x) readLines(file.path(dir, paste0(x, ".arms.tsv")))
  expect_identical(rd("a"), rd("b"))
  a <- read.delim(file.path(dir, "a.arms.tsv")); c2 <- read.delim(file.path(dir, "c.arms.tsv"))
  # the seed also drives each tumour's density fit, so the scores move a little; the p-values move more
  expect_equal(a$score_gain, c2$score_gain, tolerance = 0.01)
  expect_false(identical(a$p_gain, c2$p_gain))
})

test_that("arm_test parameters are validated; custom build warns", {
  d <- armtest_toy(armtest_passengers(6), seed = 3)
  expect_error(armtest_run(d, arm_n_perm = 0L), "arm_n_perm")
  expect_error(armtest_run(d, arm_min_markers = 0L), "arm_min_markers")
  skip_if_not_installed("jsonlite")
  f <- rae(d, build = "custom", params = fast(arm_test = TRUE, arm_n_perm = 100L), out = tempfile("armc_"))
  js <- jsonlite::fromJSON(paste(f$manifest, collapse = "\n"))
  expect_true("arm_test: no centromere coordinates for build 'custom'; each chromosome is one arm" %in% js$warnings)
  expect_identical(f$arms$arm_label, as.character(1:6))
  # with focal_only as well: one warning naming both switches
  b <- rae(d, build = "custom", params = fast(arm_test = TRUE, arm_n_perm = 100L, focal_only = TRUE), out = tempfile("armc_"))
  w <- jsonlite::fromJSON(paste(b$manifest, collapse = "\n"))$warnings
  expect_identical(grep("each chromosome is one arm", w, value = TRUE),
                   "focal_only, arm_test: no centromere coordinates for build 'custom'; each chromosome is one arm")
})

test_that("CLI and R package agree with --arm-test --arm-perm 500", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  set.seed(12)
  alt <- armtest_passengers(20, avoid = "2p")
  for (j in 1:8) alt[[j]] <- c(alt[[j]], `2p` = -0.6)
  dir <- tempfile("armapp_"); dir.create(dir)
  sf <- file.path(dir, "c.seg")
  write.table(format(armtest_toy(alt, seed = 12), scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t",
              quote = FALSE, row.names = FALSE)
  rae(read.delim(sf), build = "hg38", params = fast(arm_test = TRUE, arm_n_perm = 500L, arm_min_markers = 20L),
      out = file.path(dir, "r.rae"))
  system2(cli, c(sf, "--build", "hg38", "--n-perm", "400", "--recleave-every", "100", "--arm-test", "--arm-perm", "500",
                 "--arm-min-markers", "20", "--out", file.path(dir, "c.rae")), stdout = FALSE, stderr = FALSE)
  for (t in c("regions", "roi", "samples", "segments", "calls", "null", "arms"))
    expect_identical(readLines(file.path(dir, paste0("c.rae.", t, ".tsv"))),
                     readLines(file.path(dir, paste0("r.rae.", t, ".tsv"))), info = t)
})

test_that("the arm-state hook returns the per-tumour arm scores the arm test averages", {
  set.seed(13)
  alt <- armtest_passengers(12, avoid = "2p")
  for (j in 1:6) alt[[j]] <- c(alt[[j]], `2p` = -0.6)
  d <- armtest_toy(alt, seed = 13)
  g <- armtest_genome()
  bins <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    s <- seq(0, g$len[i] - 1e6, by = 1e6); data.frame(chrom = as.character(g$chrom[i]), start = s + 1, end = s + 1e6) }))
  p <- fast(arm_test = TRUE, arm_n_perm = 100L, arm_min_markers = 20L)
  a <- rae(d, bins = bins, build = "hg38", params = p, out = tempfile("armst_"))$arms
  st <- rae3:::rae_arm_states_cpp(d, bins, NULL, "hg38", p)
  expect_equal(nrow(st), 12 * 12)
  expect_identical(names(st), c("ID", "chrom", "arm", "eligible", "level", "coverage", "A0", "A1", "D0", "D1"))
  el <- st[st$eligible, ]
  lab <- paste0(el$chrom, el$arm)
  expect_equal(as.numeric(tapply(el$A0, lab, mean)[a$arm_label]), a$A0, tolerance = 1e-5)
  expect_equal(as.numeric(tapply(el$D0 >= 0.9, lab, mean)[a$arm_label]), a$frac_loss, tolerance = 1e-5)
  expect_equal(sum(el$D0[lab == "2p"] >= 0.9), 6)
  expect_error(rae3:::rae_arm_states_cpp(d, bins, NULL, "hg38", fast()), "arm_test or focal_only")
})
