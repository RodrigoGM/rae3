# Unit tests of single steps against an independent pure-R statement of the paper's equations.

test_that("soft discriminators follow Eq. 2", {
  a <- seq(-2, 3, by = 0.07)
  E <- 0.16; b <- 58.91
  expect_equal(rae_sigmoid(a, E, b), 1 / (exp(-b * (a - E)) + 1))
  expect_equal(rae_sigmoid(a, -0.127, -76.2, "D0"), 1 / (exp(76.2 * (a + 0.127)) + 1))
  ref <- ifelse(a > 0.3, 2 / (exp(-(a - 0.3) / log(2)) + 1) - 1, 0)
  expect_equal(rae_sigmoid(a, 0.3, 1 / log(2), "A1"), ref)
  expect_true(all(rae_sigmoid(a, 0.3, 1 / log(2), "A1") >= 0))
})

test_that("derivative matches Methods S1-C and numerical differentiation", {
  a <- seq(-1, 1, by = 0.013); E <- 0.2; b <- 30; h <- 1e-6
  num <- (rae_sigmoid(a + h, E, b) - rae_sigmoid(a - h, E, b)) / (2 * h)
  expect_equal(rae_dsigmoid(a, E, b), abs(b) * exp(b * (a - E)) / (exp(b * (a - E)) + 1)^2)
  expect_equal(rae_dsigmoid(a, E, b), num, tolerance = 1e-6)
  expect_equal(rae_dsigmoid(a, -E, -b, "D0"), abs((rae_sigmoid(a + h, -E, -b) - rae_sigmoid(a - h, -E, -b)) / (2 * h)),
               tolerance = 1e-6)
  a1 <- c(0.5, 1, 2)
  expect_equal(rae_dsigmoid(a1, 0.3, 1.44, "A1"),
               (rae_sigmoid(a1 + h, 0.3, 1.44, "A1") - rae_sigmoid(a1 - h, 0.3, 1.44, "A1")) / (2 * h), tolerance = 1e-6)
  expect_equal(rae_dsigmoid(0.1, 0.3, 1.44, "A1"), 0)
})

test_that("density reproduces stats::density and bw.nrd0", {
  set.seed(1)
  x <- c(rnorm(5000, 0, 0.03), rnorm(1500, -0.4, 0.04), rnorm(300, 0.28, 0.03))
  expect_equal(rae_bw_nrd0(x), stats::bw.nrd0(x))
  bw <- 2 * stats::bw.nrd0(x)
  ref <- stats::density(x, bw = bw, from = -1, to = 1, n = 2048)
  got <- rae_density(x, bw, -1, 1, 2048)
  expect_equal(got$x, ref$x)
  expect_lt(max(abs(got$y - ref$y)) / max(ref$y), 2e-3)
  expect_equal(which.max(got$y), which.max(ref$y))
})

test_that("BH and empirical p-values", {
  set.seed(2)
  p <- c(runif(200), runif(20, 0, 1e-4), 0.5, 0.5)
  expect_equal(rae_bh(p), stats::p.adjust(p, "BH"))
  null <- c(0.1, 0.2, 0.2, 0.3, 0.9)
  expect_equal(rae_empirical_p(c(0.05, 0.2, 0.95, 0.9), null), c(6, 5, 1, 2) / 6)   # null scores >= 0.2: four, so (4 + 1) / 6
})

test_that("generator is deterministic and stream-separated", {
  expect_identical(rae_rng_unif(5, 1, 0), rae_rng_unif(5, 1, 0))
  expect_false(identical(rae_rng_unif(5, 1, 0), rae_rng_unif(5, 1, 1)))
  u <- rae_rng_unif(20000, 7)
  expect_true(all(u >= 0 & u < 1))
  expect_lt(abs(mean(u) - 0.5), 0.01)
})

# A tumour shaped like Fig. S5A of the paper: sharp diploid peak, a gain peak near +0.28, a small
# loss peak near -0.22 and a large one at -0.4 that lies outside 2 x DN.
# These test the noise-model placement of the text (Figure S5), so they use the taylor2008 preset; the
# default (sequencing) places the detectors on the copy-number ladder instead (test-sequencing.R).
T08 <- rae_params(preset = "taylor2008")
fig_s5a <- function() {
  list(mean = c(0, 0.01, -0.01, 0.005, 0.277, -0.22, -0.40, -0.40),
       sd = rep(0.2, 8), nmark = c(60000, 50000, 40000, 30000, 6000, 3000, 30000, 20000))
}

test_that("parameterisation: E is midway between peak and half-maximum, S(H) = 0.001", {
  s <- fig_s5a()
  p <- rae_parameterize(s$mean, s$sd, s$nmark, dn = 0.19, q01 = -0.9, q025 = -0.75, params = T08)
  L <- log(999)
  expect_equal(unname(p$E["A0"]), (p$peak_gain + p$H_gain) / 2)
  expect_equal(unname(p$E["D0"]), (p$peak_loss + p$H_loss) / 2)
  expect_equal(unname(p$beta["A0"]), L / (p$E[["A0"]] - p$H_gain))
  expect_equal(unname(p$beta["D0"]), L / (p$E[["D0"]] - p$H_loss))
  expect_equal(rae_sigmoid(p$H_gain, p$E[["A0"]], p$beta[["A0"]]), 1e-3)
  expect_equal(rae_sigmoid(p$H_loss, p$E[["D0"]], p$beta[["D0"]], "D0"), 1e-3)
  # the peak window (0.5-2 x DN) picks the small -0.22 peak, not the large -0.4 one, as in Fig. S5A
  expect_lt(abs(p$peak_loss - (-0.22)), 0.02)
  expect_lt(abs(p$peak_gain - 0.277), 0.02)
  expect_true(p$beta[["D0"]] < 0 && p$beta[["A0"]] > 0)
  # same order of magnitude as the labels printed in Fig. S5A (E 0.16 / -0.127)
  expect_lt(abs(p$E[["A0"]] - 0.16), 0.03)
  expect_lt(abs(p$E[["D0"]] - (-0.127)), 0.03)
})

test_that("parameterisation: fallback fixes beta at 20", {
  s <- fig_s5a()
  p <- rae_parameterize(s$mean[1:4], s$sd[1:4], s$nmark[1:4], dn = 0.19, q01 = -0.6, q025 = -0.5, params = T08)
  expect_equal(unname(p$beta["A0"]), 20)
  expect_equal(unname(p$beta["D0"]), -20)
  expect_equal(unname(p$E["A0"]), p$H_gain + log(999) / 20)
  expect_equal(unname(p$E["D0"]), p$H_loss - log(999) / 20)
  expect_match(p$flags, "A0_fallback")
  p15 <- rae_parameterize(s$mean[1:4], s$sd[1:4], s$nmark[1:4], dn = 0.19, q01 = -0.6, q025 = -0.5,
                          params = rae_params(preset = "taylor2008", beta_fallback_loss = 15))
  expect_equal(unname(p15$beta["D0"]), -15)   # the value printed in Fig. S5C
})

test_that("parameterisation: D1 and A1 rules", {
  s <- fig_s5a()
  p <- rae_parameterize(s$mean, s$sd, s$nmark, dn = 0.19, q01 = -0.9, q025 = -0.75, params = T08)
  expect_equal(unname(p$beta["D1"]), unname(p$beta["D0"]))
  expect_equal(unname(p$E["D1"]), -0.75 - p$shift)                    # 2.5th percentile, mode-centred
  expect_equal(unname(p$beta["A1"]), 1 / log(2))
  expect_equal(unname(p$E["A1"]), max(2 * p$E[["A0"]] - p$H_gain, 1.4 * p$E[["A0"]]))
  # 2.5th percentile not negative enough -> 1st percentile
  p2 <- rae_parameterize(s$mean, s$sd, s$nmark, dn = 0.19, q01 = -0.9, q025 = -0.05, params = T08)
  expect_equal(unname(p2$E["D1"]), -0.9 - p2$shift)
  # neither is -> 1.4 x E_D0 (Table S5: E_D1 / E_D0 = 1.40 in the extreme rows)
  p3 <- rae_parameterize(s$mean, s$sd, s$nmark, dn = 0.19, q01 = -0.06, q025 = -0.05, params = T08)
  expect_equal(unname(p3$E["D1"]), 1.4 * p3$E[["D0"]])
  expect_match(p3$flags, "D1_bound")
})

test_that("parameterisation: mode centring removes a global offset", {
  s <- fig_s5a()
  p0 <- rae_parameterize(s$mean, s$sd, s$nmark, dn = 0.19, q01 = -0.9, q025 = -0.75, params = T08)
  p1 <- rae_parameterize(s$mean + 0.13, s$sd, s$nmark, dn = 0.19, q01 = -0.9 + 0.13, q025 = -0.75 + 0.13, params = T08)
  expect_equal(p1$shift - p0$shift, 0.13, tolerance = 2e-3)
  expect_equal(p1$E, p0$E, tolerance = 5e-3)
})

test_that("unified breakpoint profile rules", {
  # 20 markers; sample 1 breaks at 6 and 11; sample 2 at 6 and 15
  r <- rae_ubp(list(c(1, 6, 11), c(1, 6, 15)), list(c(5, 5, 10), c(5, 9, 6)), 20)
  expect_equal(r$start, c(1, 6, 11, 15)); expect_equal(r$end, c(5, 10, 14, 20))
  # a segment of fewer than 3 markers does not contribute its start
  r <- rae_ubp(list(c(1, 6, 8)), list(c(5, 2, 13)), 20)
  expect_equal(r$start, c(1, 8))
  # isolated single-marker region merges into the 3' region: starts 1, 10, 11 -> [10,10] joins [11,20]
  r <- rae_ubp(list(c(1, 10), c(1, 11)), list(c(9, 11), c(10, 10)), 20)
  expect_equal(r$start, c(1, 10)); expect_equal(r$end, c(9, 20))
  # contiguous single-marker regions are consolidated as one independent region: starts 10, 11, 12
  r <- rae_ubp(list(c(1, 10), c(1, 11), c(1, 12)), list(c(9, 11), c(10, 10), c(11, 9)), 20)
  expect_equal(r$start, c(1, 10, 12)); expect_equal(r$end, c(9, 11, 20))
})

test_that("regions of interest: events, peaks, plateaus, error-bar merging", {
  P <- rae_params()
  # flat significant plateau (arm-level event): one event, no peak
  r <- rae_roi(rep(0.3, 6), rep(0.01, 6), rep(1e-5, 6))
  expect_equal(nrow(r), 1); expect_false(r$is_peak)
  # singleton
  r <- rae_roi(c(0.01, 0.4, 0.01), rep(0.01, 3), c(1, 1e-5, 1))
  expect_equal(nrow(r), 1); expect_equal(c(r$first, r$last), c(2, 2))
  # one clear peak inside a broad event; neighbour within the peak's error bar is merged
  sc <- c(0.30, 0.31, 0.50, 0.495, 0.31, 0.30)
  r <- rae_roi(sc, rep(0.01, 6), rep(1e-5, 6))
  pk <- r[r$is_peak, ]
  expect_equal(nrow(pk), 1); expect_equal(c(pk$first, pk$last, pk$apex), c(3, 4, 3))
  # two peaks separated by a valley deeper than the shoulder (2 x median error)
  sc <- c(0.30, 0.50, 0.30, 0.45, 0.30)
  r <- rae_roi(sc, rep(0.01, 5), rep(1e-5, 5))
  expect_equal(sum(r$is_peak), 2)
  # a bump smaller than the shoulder is not a peak
  sc <- c(0.30, 0.50, 0.47, 0.48, 0.30)
  r <- rae_roi(sc, rep(0.01, 5), rep(1e-5, 5))
  expect_equal(sum(r$is_peak), 1)
  # regions less significant than peak_q cannot be peaks; events never cross chromosomes
  r <- rae_roi(c(0.3, 0.5, 0.3), rep(0.01, 3), rep(0.005, 3))
  expect_equal(sum(r$is_peak), 0)
  r <- rae_roi(rep(0.3, 4), rep(0.01, 4), rep(1e-5, 4), chrom = c(1, 1, 2, 2))
  expect_equal(nrow(r), 2)
})

test_that("cgh preset (released code): parameter rules of RAE beta 0.91", {
  s <- fig_s5a()
  P <- rae_params(preset = "cgh")
  p <- rae_parameterize(s$mean, s$sd, s$nmark, dn = 0.19, q01 = -0.9, q025 = -0.75, params = P)
  # loss: E is midway between the peak and zero (not the half-maximum point); gain as in the paper
  expect_equal(unname(p$E["D0"]), p$peak_loss / 2)
  expect_equal(unname(p$E["A0"]), (p$peak_gain + p$H_gain) / 2)
  expect_gte(p$H_gain, 0.025); expect_lte(p$H_loss, -0.025)          # half-maximum floor
  expect_gte(abs(p$beta[["D0"]]), 10); expect_gte(abs(p$beta[["A0"]]), 10)
  # Fig. S5A prints E = -0.127 for a loss peak near -0.25: (peak + 0)/2
  p25 <- rae_parameterize(c(s$mean[1:5], -0.254, -0.4, -0.4), s$sd, s$nmark, dn = 0.19, q01 = -0.9, q025 = -0.75, params = P)
  expect_lt(abs(p25$E[["D0"]] - (-0.127)), 0.01)
  # fallbacks: 20 for gain, 15 for loss (Fig. S5C)
  f <- rae_parameterize(s$mean[1:4], s$sd[1:4], s$nmark[1:4], dn = 0.19, q01 = -0.6, q025 = -0.5, params = P)
  expect_equal(unname(f$beta["A0"]), 20); expect_equal(unname(f$beta["D0"]), -15)
  # a very clean tumour hits the floor E = 0.025 + log(999)/100 = 0.094, the minimum in Table S5
  expect_equal(unname(f$E["A0"]), 0.025 + log(999) / 20)
  g <- rae_parameterize(c(0, 0, 0, 0.05), rep(0.2, 4), c(60000, 50000, 40000, 3000), dn = 0.09, q01 = -0.5, q025 = -0.4, params = P)
  expect_equal(unname(g$E["A0"]), 0.025 + log(999) / 100, tolerance = 0.12)
  expect_gte(g$E[["A0"]], 0.094 - 1e-9)
})

test_that("cgh preset (released code): coded peak detector merges within a two-sided doubled error", {
  P <- rae_params(preset = "cgh")
  sc <- c(0.30, 0.31, 0.50, 0.495, 0.31, 0.30)
  r <- rae_roi(sc, rep(0.01, 6), rep(1e-5, 6), params = P)
  pk <- r[r$is_peak, ]
  expect_equal(nrow(pk), 1); expect_equal(c(pk$first, pk$last, pk$apex), c(3, 4, 3))
  # flat plateau: the trace never rises above the floor, so no peak
  r <- rae_roi(rep(0.3, 6), rep(0.01, 6), rep(1e-5, 6), params = P)
  expect_equal(sum(r$is_peak), 0)
  expect_error(rae_params(preset = "nonsense"))
})
