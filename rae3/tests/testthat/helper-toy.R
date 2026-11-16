# Shared synthetic cohort for the end-to-end and plotting tests.

toy_cohort <- function(n = 40, seed = 3, binsize = 1e5, nbins = 400, chroms = 1:4) {
  set.seed(seed)
  rows <- list()
  for (j in seq_len(n)) {
    id <- sprintf("S%02d", j)
    for (ch in chroms) {
      brk <- c(0, nbins)
      val <- 0
      if (ch == 2 && j <= 16) { brk <- c(0, 150, 170, nbins); val <- c(0, -0.7, 0) }   # recurrent focal loss
      if (ch == 3 && j %% 5 == 0) { brk <- c(0, 200, 215, nbins); val <- c(0, 1.5, 0) } # focal amplification, 20 %
      if (ch == 4 && j %% 4 == 0) {                                                     # non-recurrent passengers
        a <- sample(20:300, 1); brk <- c(0, a, a + 40, nbins); val <- c(0, sample(c(-0.6, 0.5), 1), 0)
      }
      for (k in seq_along(val))
        rows[[length(rows) + 1]] <- data.frame(ID = id, chrom = ch, loc.start = brk[k] * binsize + 1,
          loc.end = brk[k + 1] * binsize, num.mark = brk[k + 1] - brk[k],
          seg.mean = val[k] + rnorm(1, 0, 0.005), seg.sd = 0.12)
    }
  }
  do.call(rbind, rows)
}
fast <- function(...) rae_params(n_perm = 400L, recleave_every = 100L, ...)
