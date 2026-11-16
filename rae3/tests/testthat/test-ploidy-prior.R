# Soft ploidy prior for the "sequencing" ladder (gap register G-14a, G-14b). Opt-in; off changes nothing.

# Samples whose segment copy numbers are exact counts (cn_counts: named by copy number, per 20 segments),
# so segment means sit on the ladder r(c) = (p c + 2(1-p)) / (p psi + 2(1-p)), mean-normalised, with only
# small noise. Purity 0.45 leaves room for the ladder to start more than one step above zero (G-13c).
ladder_counts <- function(cn_counts, n = 6, purity = 0.45, seed = 5, binsize = 1e5) {
  set.seed(seed)
  cn_one <- rep(as.numeric(names(cn_counts)), cn_counts)
  psi <- mean(cn_one)                       # 20 equal-length segments
  rows <- list()
  for (j in seq_len(n)) {
    cn <- sample(cn_one)
    chrom <- rep(1:4, each = 5); k <- rep(0:4, 4)
    r <- (purity * cn + 2 * (1 - purity)) / (purity * psi + 2 * (1 - purity))
    rows[[j]] <- data.frame(ID = sprintf("L%02d", j), chrom = chrom, loc.start = k * 80 * binsize + 1,
                            loc.end = (k + 1) * 80 * binsize, num.mark = 80,
                            seg.mean = log2(r) + rnorm(20, 0, 0.003), seg.sd = 0.1)
  }
  list(seg = do.call(rbind, rows), psi = psi,
       truth = data.frame(ID = sprintf("L%02d", seq_len(n)), purity = purity, ploidy = psi))
}
run_prior <- function(d, ...) rae(d$seg, build = "custom", params = fast(...))

test_that("prior off is identical to no prior", {
  a <- rae(toy_cohort(n = 20), build = "custom", params = fast())
  b <- rae(toy_cohort(n = 20), build = "custom", params = fast(ploidy_prior = FALSE))
  for (t in c("samples", "regions", "roi", "calls", "segments")) expect_identical(a[[t]], b[[t]], info = t)
  expect_false(any(c("ploidy_hard", "ploidy_alt", "ploidy_prior_ratio") %in% names(a$samples)))
})

test_that("prior can only raise the base copy number", {
  res <- rae(toy_cohort(), build = "custom", params = fast(ploidy_prior = TRUE))    # no resolvable ladder: columns, NA
  expect_equal(tail(names(res$samples), 4), c("flags", "ploidy_hard", "ploidy_alt", "ploidy_prior_ratio"))
  for (cn in list(c("2" = 15, "3" = 3, "4" = 2), c("2" = 5, "3" = 5, "4" = 4, "5" = 4, "6" = 2), c("1" = 5, "2" = 6, "3" = 5, "4" = 4))) {
    s <- run_prior(ladder_counts(cn), ploidy_prior = TRUE)$samples
    dp <- s$ploidy - s$ploidy_hard
    expect_true(all(dp > -1e-4 & abs(dp - round(dp)) < 1e-4))
  }
})

test_that("with the defaults a sample whose hard ploidy is about 1.3 is raised by one copy", {
  d <- ladder_counts(c("2" = 15, "3" = 3, "4" = 2))      # true 2.35; lowest level taken as one copy: 1.35
  s <- run_prior(d, ploidy_prior = TRUE)$samples
  expect_true(all(s$copy_step_source == "comb"))
  expect_equal(s$ploidy_hard, rep(d$psi - 1, nrow(s)), tolerance = 0.1)
  expect_equal(s$ploidy - s$ploidy_hard, rep(1, nrow(s)), tolerance = 1e-4)
  expect_equal(s$ploidy, rep(d$psi, nrow(s)), tolerance = 0.1)
  # the runner-up is a different candidate (here two copies up, next to the doubled component's mean)
  expect_true(all(s$ploidy_alt > s$ploidy & s$ploidy_prior_ratio > 1))
})

test_that("a high whole-genome-doubling weight raises a doubled genome with no copy-1 level", {
  d <- ladder_counts(c("2" = 5, "3" = 5, "4" = 4, "5" = 4, "6" = 2))   # true 3.65; hard 2.65
  s <- run_prior(d, ploidy_prior = TRUE, ploidy_prior_wgd = 0.8)$samples
  expect_equal(s$ploidy_hard, rep(2.65, nrow(s)), tolerance = 0.3)
  expect_equal(s$ploidy, rep(3.6, nrow(s)), tolerance = 0.3)
})

test_that("with the defaults a sample whose hard ploidy is about 2.4 is not raised", {
  d <- ladder_counts(c("1" = 5, "2" = 6, "3" = 5, "4" = 4))             # true = hard = 2.4
  s <- run_prior(d, ploidy_prior = TRUE)$samples
  expect_equal(s$ploidy_hard, rep(2.4, nrow(s)), tolerance = 0.1)
  expect_equal(s$ploidy, s$ploidy_hard, tolerance = 1e-4)
})

test_that("supplied purity/ploidy ignores the prior", {
  d <- ladder_counts(c("2" = 15, "3" = 3, "4" = 2))
  s <- rae(d$seg, build = "custom", params = fast(ploidy_prior = TRUE), purity_ploidy = d$truth)$samples
  expect_true(all(s$copy_step_source == "input"))
  expect_true(all(is.na(s$ploidy_hard) & is.na(s$ploidy_alt) & is.na(s$ploidy_prior_ratio)))
  for (n in c("ploidy_hard", "ploidy_alt", "ploidy_prior_ratio")) expect_type(s[[n]], "double")   # all NA, numeric
  expect_false(any(grepl("ploidy_ambiguous", s$flags)))
  expect_equal(s$ploidy, rep(d$psi, nrow(s)), tolerance = 1e-4)
})

test_that("ambiguous tie-breaks are flagged", {
  d <- ladder_counts(c("2" = 5, "3" = 5, "4" = 4, "5" = 4, "6" = 2))   # defaults: 2.65 vs 3.65 is a near tie
  s <- run_prior(d, ploidy_prior = TRUE)$samples
  expect_true(all(s$ploidy_prior_ratio < 3))
  expect_true(all(grepl("ploidy_ambiguous", s$flags)))
  s2 <- run_prior(d, ploidy_prior = TRUE, ploidy_prior_ambig = 1)$samples
  expect_false(any(grepl("ploidy_ambiguous", s2$flags)))
})

test_that("run.json counts the ambiguous samples only when the prior is on", {
  skip_if_not_installed("jsonlite")
  d <- ladder_counts(c("2" = 5, "3" = 5, "4" = 4, "5" = 4, "6" = 2), n = 8)
  fb <- function(res) jsonlite::fromJSON(paste(res$manifest, collapse = "\n"))$inferred_fallbacks_fired
  on <- run_prior(d, ploidy_prior = TRUE)
  expect_equal(fb(on)$`G-14b_ploidy_ambiguous_samples`, sum(grepl("ploidy_ambiguous", on$samples$flags)))
  expect_equal(fb(on)$`G-14b_ploidy_ambiguous_samples`, 8)
  expect_equal(fb(run_prior(d, ploidy_prior = TRUE, ploidy_prior_ambig = 1))$`G-14b_ploidy_ambiguous_samples`, 0)
  expect_false("G-14b_ploidy_ambiguous_samples" %in% names(fb(run_prior(d))))
})

test_that("invalid prior parameters stop the run; underflowing densities give NA, never inf", {
  d <- ladder_counts(c("2" = 15, "3" = 3, "4" = 2))
  expect_error(run_prior(d, ploidy_prior = TRUE, ploidy_prior_sd = 0), "ploidy prior needs")
  expect_error(run_prior(d, ploidy_prior = TRUE, ploidy_prior_wgd = 1.5), "ploidy prior needs")
  expect_error(run_prior(d, ploidy_prior = TRUE, ploidy_prior_ambig = 0.5), "ploidy prior needs")
  expect_error(run_prior(d, ploidy_prior = TRUE, ploidy_prior_mu1 = NaN), "finite ploidy_prior_mu1")
  expect_error(run_prior(d, ploidy_prior = TRUE, ploidy_prior_mu2 = Inf), "finite ploidy_prior_mu1")
  expect_no_error(run_prior(d, ploidy_prior = FALSE, ploidy_prior_mu1 = NaN))
  expect_no_error(run_prior(d, ploidy_prior = FALSE, ploidy_prior_sd = 0))     # only checked when on
  s <- run_prior(d, ploidy_prior = TRUE, ploidy_prior_sd = 0.005)$samples
  expect_false(any(is.infinite(s$ploidy_prior_ratio)))
  expect_true(all(is.na(s$ploidy_prior_ratio)))
})

test_that("the ratio is never infinite across small SDs; a prior with no information is flagged", {
  d <- ladder_counts(c("2" = 15, "3" = 3, "4" = 2))
  # the SD range in which the runner-up density passes through denormal values (not guaranteed to hit one)
  for (sd in c(0.01, 0.015, 0.02, 0.025, 0.03, 0.05)) {
    s <- run_prior(d, ploidy_prior = TRUE, ploidy_prior_sd = sd)$samples
    expect_false(any(is.infinite(s$ploidy_prior_ratio)), info = sd)
  }
  # every density underflows to 0: the hard choice stands, ratio NA, flagged
  s <- run_prior(d, ploidy_prior = TRUE, ploidy_prior_sd = 1e-4)$samples
  expect_equal(s$ploidy, s$ploidy_hard, tolerance = 1e-4)
  expect_true(all(is.na(s$ploidy_prior_ratio)))
  expect_true(all(grepl("ploidy_ambiguous", s$flags)))
  # a decisive underflow of the runner-up only (best density positive) is NA but not flagged
  s2 <- run_prior(d, ploidy_prior = TRUE, ploidy_prior_sd = 0.05)$samples
  expect_false(any(grepl("ploidy_ambiguous", s2$flags)))
})

test_that("rae_read() reads a run saved without the new parameters", {
  dir <- tempfile("pp_"); dir.create(dir)
  rae(toy_cohort(n = 12), build = "custom", params = fast(), out = file.path(dir, "x"))
  js <- file.path(dir, "x.run.json")
  txt <- gsub(', "ploidy_prior[a-z0-9_]*": [^,}]+', "", readLines(js))
  expect_false(any(grepl("ploidy_prior", txt)))
  writeLines(txt, js)
  r <- rae_read(file.path(dir, "x"))
  expect_false(r$params$ploidy_prior)
})

test_that("the command-line tool matches with --ploidy-prior --ploidy-prior-sd 0.3", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  d <- ladder_counts(c("2" = 15, "3" = 3, "4" = 2), n = 10)
  dir <- tempfile("rapp_"); dir.create(dir)
  sf <- file.path(dir, "c.seg")
  write.table(format(d$seg, scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t", quote = FALSE, row.names = FALSE)
  rae(read.delim(sf), build = "custom", params = fast(ploidy_prior = TRUE, ploidy_prior_sd = 0.3, ploidy_prior_ambig = 2), out = file.path(dir, "r.rae"))
  system2(cli, c(sf, "--build", "custom", "--n-perm", "400", "--recleave-every", "100", "--ploidy-prior",
                 "--ploidy-prior-sd", "0.3", "--ploidy-prior-ambig", "2", "--out", file.path(dir, "c.rae")), stdout = FALSE, stderr = FALSE)
  for (t in c("regions", "roi", "samples", "segments", "calls"))
    expect_identical(readLines(file.path(dir, paste0("c.rae.", t, ".tsv"))),
                     readLines(file.path(dir, paste0("r.rae.", t, ".tsv"))), info = t)
  expect_true("ploidy_hard" %in% names(read.delim(file.path(dir, "c.rae.samples.tsv"))))
})
