# Cohorts for the peak-independence tests (test-independence.R) and the cross-feature test (test-combo.R).

# hg38 toy cohort: chr1, 2, 4, 5 (60 Mb each) neutral apart from non-recurrent passengers; toy chr3 (60-120 Mb,
# across the centromere). `gain` samples sit at +0.3 over the whole of toy chr3. Each element of `ev` is a list
# (s = samples, a, b = bp, v = level added on [a, b]); overlapping events add up.
ind_toy <- function(n = 40, seed = 5, bs = 2.5e5, gain = integer(0), ev = list()) {
  set.seed(seed)
  rows <- list()
  add <- function(id, ch, a, b, v) rows[[length(rows) + 1]] <<- data.frame(ID = id, chrom = ch,
    loc.start = a + 1, loc.end = b, num.mark = round((b - a) / bs), seg.mean = v + rnorm(1, 0, 0.005), seg.sd = 0.12)
  for (j in seq_len(n)) {
    id <- sprintf("P%02d", j)
    for (ch in c(1, 2, 4, 5)) {
      if (j %% 4 == 0 && ch == 1 + (j %% 3)) {
        a <- sample(10:40, 1) * 1e6
        add(id, ch, 0, a, 0); add(id, ch, a, a + 4e6, sample(c(-0.6, 0.5), 1)); add(id, ch, a + 4e6, 60e6, 0)
      } else add(id, ch, 0, 60e6, 0)
    }
    v <- if (j %in% gain) 0.3 else 0
    mine <- Filter(function(e) j %in% e$s, ev)
    brk <- sort(unique(c(60e6, 120e6, unlist(lapply(mine, function(e) c(e$a, e$b))))))
    for (k in seq_len(length(brk) - 1))
      add(id, 3, brk[k], brk[k + 1], v + sum(vapply(mine, function(e) if (e$a <= brk[k] && e$b >= brk[k + 1]) e$v else 0, 0)))
  }
  do.call(rbind, rows)
}
# one amplicon at 100-110 Mb in samples 1-12, highest at 101-103 Mb, with a lower maximum at 106-108 Mb
shoulder_toy <- function() ind_toy(ev = list(list(s = 1:12, a = 100e6, b = 110e6, v = 0.6),
                                             list(s = 1:12, a = 101e6, b = 103e6, v = 0.8),
                                             list(s = 1:12, a = 106e6, b = 108e6, v = 0.5)))
# the liposarcoma pattern: a shared gain of toy chr3 (samples 1-24) carrying two amplicons 5 Mb apart in disjoint
# samples (100-102 Mb in 1-8, 107-109 Mb in 9-14), joined into one event by a lower gain between them (15-20)
lipo_toy <- function() ind_toy(gain = 1:24, ev = list(list(s = 1:8, a = 100e6, b = 102e6, v = 1.2),
                                                      list(s = 9:14, a = 107e6, b = 109e6, v = 1.2),
                                                      list(s = 15:20, a = 102e6, b = 107e6, v = 0.6)))
