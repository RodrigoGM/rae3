# Focal-only analysis (gap register G-14d): arm-level change is removed from the per-segment scores, and the
# null model is rebuilt from the same focal scores. Opt-in; off changes nothing (test-defaults.R).

# hg38 toy cohort. chr3 is the test chromosome and spans its centromere (about 90.5-93.7 Mb); chr1, 2, 4 and 5
# (60 Mb each) stay neutral apart from non-recurrent passengers. `gain` samples sit at +0.3 over the whole
# toy chr3; `amp` samples also carry an amplification 1.2 above that level at `amp_at`; `sub` samples carry
# a further +0.3 at `sub_at`.
focal_toy <- function(n = 40, seed = 7, bs = 2.5e5, chr3 = c(60e6, 120e6), gain = 1:20, amp = 1:8,
                      amp_at = c(100e6, 102e6), sub = integer(0), sub_at = c(30e6, 60e6)) {
  set.seed(seed)
  rows <- list()
  add <- function(id, ch, a, b, v) rows[[length(rows) + 1]] <<- data.frame(ID = id, chrom = ch,
    loc.start = a + 1, loc.end = b, num.mark = round((b - a) / bs), seg.mean = v + rnorm(1, 0, 0.005), seg.sd = 0.12)
  for (j in seq_len(n)) {
    id <- sprintf("F%02d", j)
    for (ch in c(1, 2, 4, 5)) {
      if (j %% 4 == 0 && ch == 1 + (j %% 3)) {
        a <- sample(10:40, 1) * 1e6
        add(id, ch, 0, a, 0); add(id, ch, a, a + 4e6, sample(c(-0.6, 0.5), 1)); add(id, ch, a + 4e6, 60e6, 0)
      } else add(id, ch, 0, 60e6, 0)
    }
    v <- if (j %in% gain) 0.3 else 0
    ev <- list()
    if (j %in% amp) ev[[length(ev) + 1]] <- c(amp_at, v + 1.2)
    if (j %in% sub) ev[[length(ev) + 1]] <- c(sub_at, v + 0.3)
    pos <- chr3[1]
    for (e in ev[order(vapply(ev, `[`, 0, 1))]) { add(id, 3, pos, e[1], v); add(id, 3, e[1], e[2], e[3]); pos <- e[2] }
    add(id, 3, pos, chr3[2], v)
  }
  do.call(rbind, rows)
}
focal_run <- function(d, build = "hg38", ...) rae(d, build = build, params = fast(...), out = tempfile("focal_"))

test_that("focal amplification on a gained arm survives, the arm gain does not", {
  d <- focal_toy()
  f <- focal_run(d, focal_only = TRUE)
  roi <- f$roi[f$roi$direction == "gain", ]
  expect_true(any(roi$chrom == 3 & roi$start <= 102e6 & roi$end >= 100e6 + 1))
  sig <- f$regions[f$regions$A_q <= 0.1, ]
  expect_gt(nrow(sig), 0)
  expect_true(all(sig$chrom == 3 & sig$start >= 100e6 + 1 & sig$end <= 102e6))
  # with the defaults the arm gain makes chr3 significant across the arm
  r <- focal_run(d)
  c3 <- r$regions[r$regions$chrom == 3, ]
  bp <- c3$end - c3$start + 1
  expect_gt(sum(bp[c3$A_q <= 0.1]) / sum(bp), 0.95)
  expect_false("focal_excluded_frac" %in% names(r$samples))
  expect_equal(f$samples$focal_excluded_frac, rep(0, 40))
  expect_type(f$samples$focal_excluded_frac, "double")   # all 0, still numeric
})

test_that("the null is rebuilt on focal scores", {
  d <- focal_toy()
  nul <- function(...) rae3:::rae_null_summary_cpp(rae3:::.seg_df(d), "hg38", fast(...))
  f <- nul(focal_only = TRUE); r <- nul()
  # the arm gains (20 of 40 samples, all of toy chr3) leave the scores, so the null's gain tail shrinks
  expect_lt(f[["A_max"]], r[["A_max"]])
  expect_lt(f[["A_q99"]], r[["A_q99"]])
})

test_that("null.tsv geometry changes with the centromere split", {
  d <- focal_toy()
  nul <- function(res) readLines(grep("null\\.tsv$", res$files, value = TRUE))
  expect_false(identical(nul(focal_run(d, focal_only = TRUE)), nul(focal_run(d))))
})

test_that("rae_fga, fga and n_altered stay absolute under focal_only", {
  d <- focal_toy(sub = 21:30, sub_at = c(65e6, 90e6))
  f <- focal_run(d, focal_only = TRUE)$samples
  r <- focal_run(d)$samples
  expect_equal(f$rae_fga, r$rae_fga, tolerance = 1e-12)
  expect_equal(f$fga, r$fga, tolerance = 1e-12)
  expect_equal(f$n_altered, r$n_altered)          # a segment cut at the centromere counts once
  # the gained samples' toy chr3 (60-120 Mb, one segment across the centromere) is altered and was cut
  gained <- sprintf("F%02d", 9:20)
  expect_true(all(r$rae_fga[r$ID %in% gained] >= 60 / 300 - 1e-9))
  # cut at the centromere: toy chr3 and the neutral chr4 and chr5 (hg38 centromeres near 50 and 48 Mb)
  expect_true(all(f$n_segments[f$ID %in% gained] == r$n_segments[r$ID %in% gained] + 3))
})

test_that("deviations longer than focal_max_bp are excluded and reported", {
  # a 30 Mb gain inside a neutral p arm (toy chr3 is 0-120 Mb)
  d <- focal_toy(chr3 = c(0, 120e6), gain = integer(0), amp = integer(0), sub = 1:20, sub_at = c(30e6, 60e6))
  f <- focal_run(d, focal_only = TRUE)
  ex <- setNames(f$samples$focal_excluded_frac, f$samples$ID)
  expect_equal(unname(ex[sprintf("F%02d", 1:20)]), rep(30 / 360, 20), tolerance = 1e-3)
  expect_equal(unname(ex[sprintf("F%02d", 21:40)]), rep(0, 20))
  inside <- function(res) { g <- res$regions; g[g$chrom == 3 & g$start > 30e6 & g$end <= 60e6, ] }
  expect_true(all(inside(f)$A_q > 0.1))
  # the same gain is kept, and significant, when the limit is above its length
  k <- focal_run(d, focal_only = TRUE, focal_max_bp = 4e7)
  expect_true(all(k$samples$focal_excluded_frac == 0))
  expect_true(all(inside(k)$A_q <= 0.1))
})

# Samples on an exact copy-number ladder (purity 0.7): chr1, 2, 4 at two copies, chr5 at three, toy chr3
# (0-120 Mb) at one copy with a 2 Mb zero-copy region at 100-102 Mb.
homdel_toy <- function(n = 6, purity = 0.7, bs = 2.5e5, seed = 11) {
  set.seed(seed)
  segs <- data.frame(chrom = c(1, 2, 4, 5, 3, 3, 3), a = c(0, 0, 0, 0, 0, 100e6, 102e6),
                     b = c(60e6, 60e6, 60e6, 60e6, 100e6, 102e6, 120e6), cn = c(2, 2, 2, 3, 1, 0, 1))
  len <- segs$b - segs$a
  psi <- sum(segs$cn * len) / sum(len)
  r <- (purity * segs$cn + 2 * (1 - purity)) / (purity * psi + 2 * (1 - purity))
  do.call(rbind, lapply(seq_len(n), function(j)
    data.frame(ID = sprintf("H%02d", j), chrom = segs$chrom, loc.start = segs$a + 1, loc.end = segs$b,
               num.mark = round(len / bs), seg.mean = log2(r) + rnorm(nrow(segs), 0, 0.003), seg.sd = 0.1)))
}

test_that("homozygous deletion on a lost arm keeps D1", {
  f <- focal_run(homdel_toy(), focal_only = TRUE)
  expect_true(all(f$samples$copy_step_source == "comb"))
  sg <- f$segments
  hd <- sg[sg$chrom == 3 & sg$loc.start == 100e6 + 1, ]
  expect_equal(nrow(hd), 6)
  expect_true(all(hd$D1 >= 0.9))
  expect_true(all(hd$D0 >= 0.9))               # one copy below its arm
  arm <- sg[sg$chrom == 3 & sg$loc.start != 100e6 + 1, ]
  expect_true(all(arm$D0 < 0.1 & arm$D1 < 0.1))  # the arm loss itself is removed
  r <- focal_run(homdel_toy())
  expect_true(all(r$segments$D0[r$segments$chrom == 3] >= 0.9))   # without focal_only the arm is a loss
})

test_that("step-fallback samples use the log2 residual", {
  b <- data.frame(chrom = c(rep("1", 248), rep("2", 100)), start = c(0:247, 0:99) * 1e6 + 1, end = c(1:248, 1:100) * 1e6)
  # chr1 gained (log2 1.5) with an amplification; chr2 lost (log2 0.5) with a near-zero-copy segment at 40-42 Mb
  seg <- data.frame(ID = "S1", chrom = c("1", "1", "1", "2", "2", "2"),
                    loc.start = c(1, 150e6 + 1, 152e6 + 1, 1, 40e6 + 1, 42e6 + 1),
                    loc.end = c(150e6, 152e6, 248e6, 40e6, 42e6, 100e6),
                    num.mark = c(150L, 2L, 96L, 40L, 2L, 58L),
                    seg.mean = c(log2(1.5), log2(3), log2(1.5), -1, -2.3, -1), seg.sd = 0.1)
  E <- c(0.3, 1.2, -0.3, -2); B <- c(20, 8, -20, -20)
  prm <- rae_params(focal_only = TRUE, arm_baseline = "median")
  hook <- function(fallback) rae3:::rae_focal_scores_cpp(seg, b, "hg38", prm, step = 0.5, lad_m = 1, lad_r0 = 0,
                                                         ploidy = 2, fallback = fallback, E = E, B = B)
  o <- hook(TRUE)
  sig <- function(a, k) 1 / (exp(-B[k] * (a - E[k])) + 1)
  # log2: residual = mean - arm level (log2 1.5 on chr1, -1 on chr2); D1 at the absolute mean
  res <- o$mean - ifelse(o$chrom == 1, log2(1.5), -1)
  expect_equal(o$A0, sig(res, 1), tolerance = 1e-6)
  expect_equal(o$A1, sig(res, 2), tolerance = 1e-6)
  expect_equal(o$D0, sig(res, 3), tolerance = 1e-6)
  expect_equal(o$D1, sig(o$mean, 4), tolerance = 1e-6)
  hd <- o$mean < -2
  expect_equal(sum(hd), 1)
  expect_gt(o$D1[hd], 0.99)                  # absolute: a homozygous deletion
  expect_lt(sig(res[hd], 4), 0.01)           # at the residual it would not be
  # with a usable ladder the residual is in copy units: c' = kmode + (c - b), at x(max(c', 0))
  u <- hook(FALSE)
  x <- function(c) log2(pmax(0 + c * 0.5, 1e-6))
  cp <- (2^u$mean - 0) / 0.5; base <- ifelse(u$chrom == 1, 3, 1)
  expect_equal(u$A0, sig(x(pmax(2 + cp - base, 0)), 1), tolerance = 1e-6)
  expect_equal(u$D1, sig(u$mean, 4), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(u$A0, o$A0)))
})

test_that("the length rule and excluded fraction on one sample", {
  b <- data.frame(chrom = "1", start = 0:247 * 1e6 + 1, end = 1:248 * 1e6)
  # p arm (0-123 Mb): neutral, a 25 Mb gain, neutral; q arm: neutral with a 5 Mb loss
  seg <- data.frame(ID = "S1", chrom = "1", loc.start = c(1, 20e6 + 1, 45e6 + 1, 180e6 + 1, 185e6 + 1),
                    loc.end = c(20e6, 45e6, 180e6, 185e6, 248e6), num.mark = c(20L, 25L, 135L, 5L, 63L),
                    seg.mean = c(0, 0.6, 0, -0.6, 0), seg.sd = 0.1)
  E <- c(0.3, 1.2, -0.3, -2); B <- c(20, 8, -20, -20)
  prm <- rae_params(focal_only = TRUE, arm_baseline = "median")
  o <- rae3:::rae_focal_scores_cpp(seg, b, "hg38", prm, step = NA, lad_m = NA, lad_r0 = NA, ploidy = NA,
                                   fallback = FALSE, E = E, B = B)
  g <- o[o$mean > 0.5, ]
  expect_true(all(c(g$A0, g$A1, g$D0, g$D1) == 0))              # 25 Mb > 20 Mb
  expect_gt(o$D0[o$mean < -0.5], 0.99)                           # 5 Mb kept
  expect_equal(attr(o, "excluded_frac"), 25e6 / 248e6)
  prm$focal_max_bp <- 3e7
  o2 <- rae3:::rae_focal_scores_cpp(seg, b, "hg38", prm, step = NA, lad_m = NA, lad_r0 = NA, ploidy = NA,
                                    fallback = FALSE, E = E, B = B)
  expect_gt(o2$A0[o2$mean > 0.5], 0.99)
  expect_equal(attr(o2, "excluded_frac"), 0)
})

test_that("custom build runs with a warning", {
  skip_if_not_installed("jsonlite")
  f <- focal_run(focal_toy(), build = "custom", focal_only = TRUE)
  js <- jsonlite::fromJSON(paste(f$manifest, collapse = "\n"))
  expect_true("focal_only: no centromere coordinates for build 'custom'; each chromosome is one arm" %in% js$warnings)
  sig <- f$regions[f$regions$A_q <= 0.1, ]
  expect_true(all(sig$chrom == 3 & sig$start >= 100e6 + 1 & sig$end <= 102e6))
  g <- focal_run(focal_toy(), build = "custom")
  expect_false(any(grepl("focal_only", jsonlite::fromJSON(paste(g$manifest, collapse = "\n"))$warnings)))
})

test_that("focal_max_bp is validated and the column follows the ploidy-prior columns", {
  expect_error(focal_run(focal_toy(n = 10), focal_only = TRUE, focal_max_bp = 0), "focal_max_bp")
  f <- focal_run(focal_toy(n = 10, gain = 1:5, amp = 1:3), focal_only = TRUE, ploidy_prior = TRUE)
  expect_equal(tail(names(f$samples), 5), c("flags", "ploidy_hard", "ploidy_alt", "ploidy_prior_ratio", "focal_excluded_frac"))
})

test_that("CLI and R package agree with --focal-only --focal-max-mb 10", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  dir <- tempfile("fapp_"); dir.create(dir)
  sf <- file.path(dir, "c.seg")
  write.table(format(focal_toy(sub = 21:30, sub_at = c(65e6, 80e6)), scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t",
              quote = FALSE, row.names = FALSE)
  rae(read.delim(sf), build = "hg38", params = fast(focal_only = TRUE, focal_max_bp = 1e7), out = file.path(dir, "r.rae"))
  system2(cli, c(sf, "--build", "hg38", "--n-perm", "400", "--recleave-every", "100", "--focal-only",
                 "--focal-max-mb", "10", "--out", file.path(dir, "c.rae")), stdout = FALSE, stderr = FALSE)
  for (t in c("regions", "roi", "samples", "segments", "calls", "null"))
    expect_identical(readLines(file.path(dir, paste0("c.rae.", t, ".tsv"))),
                     readLines(file.path(dir, paste0("r.rae.", t, ".tsv"))), info = t)
  s <- read.delim(file.path(dir, "c.rae.samples.tsv"))
  expect_true(any(s$focal_excluded_frac > 0))
})
