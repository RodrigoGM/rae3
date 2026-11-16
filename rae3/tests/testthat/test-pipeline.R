# End-to-end tests on a small synthetic cohort.

test_that("pipeline finds the recurrent events and nothing on the null chromosomes", {
  seg <- toy_cohort()
  res <- rae(seg, build = "custom", params = fast())
  expect_s3_class(res, "rae_result")
  ev <- res$roi[res$roi$type == "region", ]
  loss2 <- ev[ev$direction == "loss" & ev$chrom == 2, ]
  expect_equal(nrow(loss2), 1)
  expect_lte(loss2$start, 15000001); expect_gte(loss2$end, 17000000)
  expect_lt(loss2$end - loss2$start, 3e6)
  gain3 <- ev[ev$direction == "gain" & ev$chrom == 3, ]
  expect_equal(nrow(gain3), 1)
  expect_equal(c(gain3$start, gain3$end), c(20000001, 21500000))
  expect_equal(gain3$freq_high, 0.2)   # called as multi-copy amplification in 8 of 40
  expect_equal(sum(ev$chrom %in% c(1, 4)), 0)
  expect_equal(nrow(res$samples), 40)
  expect_true(all(file.exists(res$files)))
  expect_equal(length(res$files), 8)
  # frequency at the apex agrees with the design (16 of 40)
  expect_equal(loss2$freq, 0.4, tolerance = 0.05)
})

test_that("same seed gives identical tables; another seed changes only the stochastic parts", {
  seg <- toy_cohort()
  a <- rae(seg, build = "custom", params = fast())
  b <- rae(seg, build = "custom", params = fast())
  expect_identical(a$regions, b$regions); expect_identical(a$roi, b$roi); expect_identical(a$samples, b$samples)
  d <- rae(seg, build = "custom", params = fast(seed = 99))
  expect_false(identical(a$regions$A_p, d$regions$A_p))
  expect_equal(a$regions$start, d$regions$start)
})

test_that("seg-only input runs with the documented fallbacks and reports them", {
  seg <- toy_cohort()[, 1:6]
  res <- rae(seg, build = "custom", params = fast())
  expect_true(all(grepl("dn_from_sd", res$samples$flags)))
  expect_true(any(grepl("G-01d_virtual_marker_grid\": true", res$manifest)))
  expect_true(any(res$roi$direction == "loss" & res$roi$chrom == 2))
})

test_that("normals flag germline regions", {
  seg <- toy_cohort()
  nrm <- seg[seg$ID %in% sprintf("S%02d", 1:6), ]
  nrm$ID <- sub("S", "N", nrm$ID)
  res <- rae(seg, normals = nrm, build = "custom", params = fast())
  hit <- res$roi[res$roi$direction == "loss" & res$roi$chrom == 2 & res$roi$type == "region", ]
  expect_gte(hit$germline_normals_max, 2)
})

test_that("non-autosomal rows are dropped and counted; bad input errors", {
  seg <- toy_cohort()
  x <- rbind(seg, transform(seg[1:3, ], chrom = "X"))
  res <- rae(x, build = "custom", params = fast())
  expect_true(any(grepl("rows_non_autosomal_dropped\": 3", res$manifest)))
  expect_error(rae(transform(seg, chrom = "Y"), build = "custom", params = fast()), "no autosomal")
  expect_error(rae_params(nonsense = 1), "unknown RAE parameter")
})

test_that("QDNAseq objects are converted and analysed like their seg export", {
  skip_if_not_installed("QDNAseq"); skip_if_not_installed("Biobase")
  set.seed(5)
  nb <- 300; ns <- 12
  bins <- data.frame(chromosome = rep(c("1", "2"), each = nb), start = rep((0:(nb - 1)) * 1e5 + 1, 2),
                     end = rep((1:nb) * 1e5, 2), bases = 100, gc = 50, mappability = 100, blacklist = 0,
                     residual = 0, use = TRUE, stringsAsFactors = FALSE)
  rownames(bins) <- sprintf("%s:%d-%d", bins$chromosome, bins$start, bins$end)
  truth <- matrix(0, 2 * nb, ns, dimnames = list(rownames(bins), sprintf("Q%02d", 1:ns)))
  truth[(nb + 100):(nb + 130), 1:5] <- -0.7
  cn <- truth + rnorm(length(truth), 0, 0.1)
  cn[sample(length(cn), 0.02 * length(cn))] <- NA
  sgm <- truth; sgm[is.na(cn)] <- NA
  obj <- new("QDNAseqCopyNumbers", bins = Biobase::AnnotatedDataFrame(bins), copynumber = 2^cn,
             phenodata = Biobase::AnnotatedDataFrame(data.frame(name = colnames(cn), row.names = colnames(cn))))
  Biobase::assayDataElement(obj, "segmented") <- 2^sgm
  inp <- as_seg(obj, bin_noise = TRUE)
  expect_equal(sort(unique(inp$seg$seg.mean)), c(-0.7, 0), tolerance = 1e-9)
  expect_equal(nrow(inp$stats), ns)
  expect_equal(median(inp$stats$dn), 0.954 * 0.1, tolerance = 0.08)   # DN of Gaussian noise = 0.954 sd
  # segments skip missing bins but keep their extent
  expect_true(all(inp$seg$num.mark <= (inp$seg$loc.end - inp$seg$loc.start + 1) / 1e5))
  a <- rae(obj, build = "custom", params = fast(), bin_noise = TRUE)
  b <- rae(inp$seg, bins = inp$bins, sample_stats = inp$stats, build = "custom", params = fast())
  expect_identical(a$regions, b$regions)
  expect_true(all(grepl("^-$|ubp_excluded|fallback|D1|A1", a$samples$flags)))   # no seg-only fallbacks fired
  expect_false(any(grepl("from_sd|from_pseudo|from_breakpoints", a$samples$flags)))
  # default (bin_noise = FALSE): the same as a plain seg, which is all the command-line tool sees
  plain <- as_seg(obj)
  expect_null(plain$stats)
  expect_false("seg.sd" %in% names(plain$seg))
  d <- rae(obj, build = "custom", params = fast())
  e <- rae(plain$seg, bins = plain$bins, build = "custom", params = fast())
  expect_identical(d$regions, e$regions)
  expect_true(any(grepl("from_breakpoints|from_pseudo|from_sd", d$samples$flags)))
})

test_that("command-line tool and R package write identical tables", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  seg <- toy_cohort()
  d <- tempfile("raecli_"); dir.create(d)
  f <- file.path(d, "toy.seg")
  write.table(format(seg, scientific = FALSE, trim = TRUE, digits = 15), f, sep = "\t", quote = FALSE, row.names = FALSE)
  seg2 <- read.delim(f)
  res <- rae(seg2, build = "custom", params = fast(), out = file.path(d, "r.rae"))
  system2(cli, c(f, "--build", "custom", "--n-perm", "400", "--recleave-every", "100", "--out", file.path(d, "c.rae")),
          stdout = FALSE, stderr = FALSE)
  for (s in c(".regions.tsv", ".roi.tsv", ".samples.tsv", ".segments.tsv", ".calls.tsv", ".null.tsv"))
    expect_identical(readLines(file.path(d, paste0("r.rae", s))), readLines(file.path(d, paste0("c.rae", s))), info = s)
})

test_that("cgh preset runs end to end and records itself in the manifest", {
  seg <- toy_cohort()
  P <- rae_params(preset = "cgh", n_perm = 100L, recleave_every = 10L)
  res <- rae(seg, build = "custom", params = P)
  expect_true(any(grepl("\"preset\": \"cgh\"", res$manifest)))
  expect_true(any(grepl("\"d0_midpoint_zero\": 1", res$manifest)))
  ev <- res$roi[res$roi$type == "region", ]
  expect_true(any(ev$direction == "loss" & ev$chrom == 2))
  expect_equal(sum(ev$chrom %in% c(1, 4)), 0)
})

test_that("population CNV overlap is reported for regions of interest and never filters", {
  seg <- toy_cohort()
  f <- tempfile(fileext = ".tsv")
  writeLines(c("chrom\tstart\tend\ttype\taf", "2\t15000001\t16000000\tloss\t0.05", "2\t15500001\t17000000\tloss\t0.02",
               "3\t1\t1000\tgain\t0.5"), f)
  a <- rae(seg, build = "custom", params = fast(), cnv = "none")
  b <- rae(seg, build = "custom", params = fast(), cnv = f)
  expect_equal(nrow(a$roi), nrow(b$roi))
  expect_true(all(is.na(a$roi$cnv_n)))
  hit <- b$roi[b$roi$direction == "loss" & b$roi$chrom == 2 & b$roi$type == "region", ]
  expect_equal(hit$cnv_n, 2); expect_equal(hit$cnv_loss_pct, 100); expect_equal(hit$likely_polymorphic, 1)
  g <- b$roi[b$roi$direction == "gain" & b$roi$chrom == 3, ]
  expect_true(all(g$cnv_n == 0)); expect_true(all(g$likely_polymorphic == 0))
})
