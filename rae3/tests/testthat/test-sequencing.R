# Preset "sequencing" (gap register G-13): detectors on the copy-number ladder, optional purity/ploidy.

# A cohort with a known absolute scale: integer copy numbers per segment, log2 ratio from purity p and the
# tumour's own mean ploidy, plus small segment-mean noise.
ladder_cohort <- function(n = 30, seed = 11, purity = 0.7, nbins = 400, binsize = 1e5) {
  set.seed(seed)
  rows <- list(); truth <- list()
  for (j in seq_len(n)) {
    id <- sprintf("T%02d", j)
    segs <- list()
    for (ch in 1:4) {
      cut <- sort(sample(20:(nbins - 20), 4)); brk <- c(0, cut, nbins)
      cn <- sample(c(1, 2, 2, 2, 3, 3, 4), length(brk) - 1, replace = TRUE)
      if (ch == 3 && j %% 3 == 0) cn[2] <- 10                      # amplified under both rules (>= 5, >= 9)
      segs[[ch]] <- data.frame(chrom = ch, s = brk[-length(brk)], e = brk[-1], cn = cn)
    }
    sg <- do.call(rbind, segs)
    psi <- sum(sg$cn * (sg$e - sg$s)) / sum(sg$e - sg$s)
    r <- (purity * sg$cn + 2 * (1 - purity)) / (purity * psi + 2 * (1 - purity))
    rows[[j]] <- data.frame(ID = id, chrom = sg$chrom, loc.start = sg$s * binsize + 1, loc.end = sg$e * binsize,
                            num.mark = sg$e - sg$s, seg.mean = log2(r) + rnorm(nrow(sg), 0, 0.01), seg.sd = 0.15)
    truth[[j]] <- data.frame(ID = id, purity = purity, ploidy = psi)
  }
  list(seg = do.call(rbind, rows), truth = do.call(rbind, truth))
}

test_that("sequencing is the default preset and differs from taylor2008 only in the ladder", {
  p <- rae_params(); q <- rae_params(preset = "taylor2008")
  expect_identical(p, rae_params(preset = "sequencing"))
  expect_equal(p$preset, "sequencing")
  expect_true(p$copy_step); expect_true(p$a1_logistic); expect_equal(p$amp_threshold, 0.5)
  diff <- names(q)[!mapply(identical, q, p[names(q)])]
  expect_setequal(diff, c("copy_step", "a1_logistic", "amp_threshold", "preset"))
  expect_false(rae_params(preset = "cgh")$copy_step)
  expect_error(rae_params(preset = "paper")); expect_error(rae_params(preset = "original"))
})

test_that("the array presets are untouched: no ladder, no new flags", {
  for (pr in c("taylor2008", "cgh")) {
    res <- rae(toy_cohort(), build = "custom", params = fast(preset = pr))
    expect_true(all(res$samples$copy_step_source == "-"), info = pr)
    expect_true(all(is.na(res$samples$purity)), info = pr)
    expect_false(any(grepl("step_fallback", res$samples$flags)), info = pr)
  }
})

test_that("rae_read() maps the preset names of older runs (paper, original)", {
  d <- tempfile("legacy_"); dir.create(d)
  res <- rae(toy_cohort(n = 12), build = "custom", params = fast(preset = "cgh"), out = file.path(d, "x"))
  js <- file.path(d, "x.run.json")
  writeLines(sub('"preset": "cgh"', '"preset": "original"', readLines(js)), js)
  expect_equal(rae_read(file.path(d, "x"))$params$preset, "cgh")
})

test_that("the comb recovers the one-copy step and places the detectors on the ladder", {
  d <- ladder_cohort()
  res <- rae(d$seg, build = "custom", params = fast(preset = "sequencing"))
  s <- merge(res$samples, d$truth, by = "ID", suffixes = c("", "_true"))
  ok <- s$copy_step_source == "comb"
  expect_gt(mean(ok), 0.8)
  step_true <- with(s, purity_true / (purity_true * ploidy_true + 2 * (1 - purity_true)))
  expect_lt(median(abs(s$copy_step[ok] / step_true[ok] - 1)), 0.1)
  expect_gt(mean(abs(s$ploidy[ok] - s$ploidy_true[ok]) < 0.5), 0.7)
  # A0 sits half a copy above the rounded ploidy, D0 half a copy below
  x <- s[ok, ][1, ]
  r <- function(c) x$purity * c + 2 * (1 - x$purity)
  k <- max(1, round(x$ploidy))
  expect_equal(x$E_A0 - x$E_D0, log2(r(k + 0.5) / r(k - 0.5)), tolerance = 1e-6)
  expect_false(any(grepl("A0_fallback|D0_fallback", s$flags[ok])))
})

test_that("supplied purity and ploidy set the ladder exactly and amplifications are called", {
  d <- ladder_cohort()
  res <- rae(d$seg, build = "custom", params = fast(preset = "sequencing"), purity_ploidy = d$truth)
  s <- merge(res$samples, d$truth, by = "ID", suffixes = c("", "_true"))
  expect_true(all(s$copy_step_source == "input"))
  expect_equal(s$purity, s$purity_true, tolerance = 1e-5)   # samples.tsv keeps 6 significant digits
  expect_equal(s$ploidy, s$ploidy_true, tolerance = 1e-5)
  # chromosome 3, second segment is 10 copies in every third tumour: called amplified (state 2) there
  cl <- res$calls
  amp_ids <- sprintf("T%02d", which(seq_len(30) %% 3 == 0))
  reg <- cl[cl$chrom == 3, ]
  amp_any <- vapply(amp_ids, function(i) any(reg[[i]] == 2, na.rm = TRUE), logical(1))
  expect_true(all(amp_any))
  others <- setdiff(sprintf("T%02d", 1:30), amp_ids)
  expect_false(any(vapply(others, function(i) any(cl[[i]] == 2, na.rm = TRUE), logical(1))))
})

test_that("purity/ploidy can come as sample_stats columns, and a flat profile falls back", {
  d <- ladder_cohort(n = 12)
  a <- rae(d$seg, build = "custom", params = fast(preset = "sequencing"), purity_ploidy = d$truth)
  b <- rae(d$seg, build = "custom", params = fast(preset = "sequencing"), sample_stats = d$truth)
  expect_identical(a$samples, b$samples)
  flat <- d$seg; flat$seg.mean[flat$ID == "T01"] <- rnorm(sum(flat$ID == "T01"), 0, 0.002)
  f <- rae(flat, build = "custom", params = fast(preset = "sequencing"))
  expect_match(f$samples$flags[f$samples$ID == "T01"], "step_fallback")
  expect_match(paste(f$manifest, collapse = ""), "G-13b_copy_step_fallback_samples")
})

test_that("the command-line tool gives the same result with --preset sequencing --purity-ploidy", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  d <- ladder_cohort(n = 12)
  dir <- tempfile("raeseq_"); dir.create(dir)
  sf <- file.path(dir, "c.seg"); pf <- file.path(dir, "pp.tsv")
  write.table(format(d$seg, scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(format(d$truth, digits = 15), pf, sep = "\t", quote = FALSE, row.names = FALSE)
  res <- rae(read.delim(sf), build = "custom", params = fast(preset = "sequencing"),
             purity_ploidy = read.delim(pf), out = file.path(dir, "r.rae"))
  # --preset after the overrides: the CLI applies it first, so --n-perm 400 still holds
  system2(cli, c(sf, "--build", "custom", "--n-perm", "400", "--recleave-every", "100", "--preset", "sequencing",
                 "--purity-ploidy", pf, "--out", file.path(dir, "c.rae")), stdout = FALSE, stderr = FALSE)
  for (t in c("regions", "roi", "samples", "segments", "calls"))
    expect_identical(readLines(file.path(dir, paste0("c.rae.", t, ".tsv"))),
                     readLines(file.path(dir, paste0("r.rae.", t, ".tsv"))), info = t)
})

test_that("the command-line tool defaults to sequencing, as rae() does", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  d <- ladder_cohort(n = 12)
  dir <- tempfile("raedef_"); dir.create(dir)
  sf <- file.path(dir, "c.seg")
  write.table(format(d$seg, scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t", quote = FALSE, row.names = FALSE)
  rae(read.delim(sf), build = "custom", params = fast(), out = file.path(dir, "r.rae"))
  system2(cli, c(sf, "--build", "custom", "--n-perm", "400", "--recleave-every", "100", "--out", file.path(dir, "c.rae")),
          stdout = FALSE, stderr = FALSE)
  for (t in c("regions", "roi", "samples", "segments", "calls"))
    expect_identical(readLines(file.path(dir, paste0("c.rae.", t, ".tsv"))),
                     readLines(file.path(dir, paste0("r.rae.", t, ".tsv"))), info = t)
  expect_match(paste(readLines(file.path(dir, "c.rae.run.json")), collapse = ""), '"preset": "sequencing"')
})

test_that("rae_plot_frequency() draws stacked and single-state frequencies", {
  skip_if_not_installed("ggplot2")
  res <- rae(toy_cohort(), build = "hg19", params = fast())
  expect_s3_class(rae_plot_frequency(res), "ggplot")
  expect_s3_class(rae_plot_frequency(res, type = "gain"), "ggplot")
  expect_s3_class(rae_plot_frequency(res, type = c("amplification", "deletion"), max_pct = 30), "ggplot")
  expect_s3_class(rae_plot_frequency(res, type = "loss", chrom = 2), "ggplot")
  expect_error(rae_plot_frequency(res, type = c("gain", "amplification")), "stacked")
  expect_error(rae_plot_frequency(res, type = c("stacked", "loss")), "combined")
  expect_error(rae_plot_frequency(res, type = "nonsense"), "type must be")
  # stacked bars reach freq_gain (top) and freq_loss (bottom) in percent
  b <- ggplot2::layer_data(rae_plot_frequency(res), 4)
  expect_lte(max(b$ymax), 100 + 0.045 * 100 + 1e-9)
})
