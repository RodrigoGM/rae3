# Fraction of the genome altered per sample: rae_fga (RAE scores), fga (|log2 ratio| > fga_threshold) and
# n_altered (gap register G-12a, G-12e).

# The tumour's segment-covered span: first segment start to last segment end on each autosome, summed.
covered_span <- function(sg) {
  sp <- tapply(seq_len(nrow(sg)), sg$chrom, function(i) max(sg$loc.end[i]) - min(sg$loc.start[i]) + 1)
  sum(as.numeric(sp))
}
seg_len <- function(sg) sg$loc.end - sg$loc.start + 1

test_that("rae_fga counts RAE-altered segments and fga counts |log2| > 0.2, both over the covered span", {
  seg <- toy_cohort(n = 20)
  seg$loc.start <- seg$loc.start + 5e5   # leave each chromosome's first 0.5 Mb unmeasured
  res <- rae(seg, build = "hg38", params = fast())
  expect_identical(names(res$samples)[names(res$samples) %in% c("rae_fga", "fga", "n_altered")],
                   c("rae_fga", "fga", "n_altered"))
  expect_false("cna_burden" %in% names(res$samples))
  sg <- res$segments
  for (id in unique(sg$ID)) {
    x <- sg[sg$ID == id, ]; len <- seg_len(x); span <- covered_span(x)
    s <- res$samples[res$samples$ID == id, ]
    by_rae <- x$A0 >= 0.75 | x$D0 >= 0.75; by_thr <- abs(x$seg.mean.centred + s$shift) > 0.2   # the ratio as supplied
    expect_equal(s$rae_fga, sum(len[by_rae]) / span, tolerance = 1e-5, info = id)
    expect_equal(s$fga, sum(len[by_thr]) / span, tolerance = 1e-5, info = id)
    expect_identical(s$n_altered, sum(by_rae), info = id)
  }
  expect_true(any(res$samples$rae_fga > 0) && any(res$samples$fga > 0))
  for (p in c("sequencing", "cgh", "taylor2008")) {
    expect_false(rae_params(preset = p)$fga_assembly, info = p)
    expect_identical(rae_params(preset = p)$fga_threshold, 0.2, info = p)
  }
})

test_that("fga_threshold moves fga only; fga does not depend on the preset", {
  seg <- toy_cohort(n = 20)
  a <- rae(seg, build = "hg38", params = fast())
  hi <- rae(seg, build = "hg38", params = fast(fga_threshold = 0.6))   # toy alterations are 0.5-1.5
  cg <- rae(seg, build = "hg38", params = fast(preset = "cgh"))
  expect_true(all(hi$samples$fga <= a$samples$fga) && any(hi$samples$fga < a$samples$fga))
  expect_identical(hi$samples$rae_fga, a$samples$rae_fga)
  expect_identical(hi$samples$n_altered, a$samples$n_altered)
  expect_identical(cg$samples$fga, a$samples$fga)
  expect_false(identical(cg$samples$shift, a$samples$shift))   # different centring, same fga
})

test_that("fga_threshold_loss sets the loss side separately; unset, losses use fga_threshold", {
  seg <- toy_cohort(n = 20)   # toy gains +0.5 and +1.5, losses -0.6 and -0.7
  a <- rae(seg, build = "hg38", params = fast(fga_threshold = 0.4))
  same <- rae(seg, build = "hg38", params = fast(fga_threshold = 0.4, fga_threshold_loss = 0.4))
  pos <- rae(seg, build = "hg38", params = fast(fga_threshold = 0.4, fga_threshold_loss = 0.65))
  neg <- rae(seg, build = "hg38", params = fast(fga_threshold = 0.4, fga_threshold_loss = -0.65))
  expect_identical(same$samples$fga, a$samples$fga)
  expect_identical(neg$samples$fga, pos$samples$fga)
  sg <- pos$segments
  for (id in unique(sg$ID)) {
    x <- sg[sg$ID == id, ]; s <- pos$samples[pos$samples$ID == id, ]
    r <- x$seg.mean.centred + s$shift
    expect_equal(s$fga, sum(seg_len(x)[r > 0.4 | r < -0.65]) / covered_span(x), tolerance = 1e-5, info = id)
  }
  expect_true(any(pos$samples$fga < a$samples$fga))      # the -0.6 losses drop out
  expect_identical(pos$samples$rae_fga, a$samples$rae_fga)
  expect_null(rae_params()$fga_threshold_loss)
})

test_that("fga_assembly divides both by the whole assembly and changes nothing else", {
  seg <- toy_cohort(n = 20)
  a <- rae(seg, build = "hg38", params = fast())
  b <- rae(seg, build = "hg38", params = fast(fga_assembly = TRUE))
  span <- covered_span(a$segments[a$segments$ID == a$samples$ID[1], ])
  for (col in c("rae_fga", "fga")) {
    k <- a$samples[[col]] > 0
    expect_true(any(k), info = col)
    r <- b$samples[[col]][k] / a$samples[[col]][k]   # every toy tumour covers the same span
    expect_equal(r, rep(r[1], sum(k)), tolerance = 1e-4, info = col)
    expect_equal(r[1], span / 2875001522, tolerance = 1e-3, info = col)   # hg38 autosomes, 2,875,001,522 bp
  }
  expect_identical(b$samples$n_altered, a$samples$n_altered)
  for (t in c("regions", "roi", "segments", "calls")) expect_identical(b[[t]], a[[t]], info = t)
})

test_that("rae_read() maps runs saved with cna_burden and no fga_assembly", {
  seg <- toy_cohort(n = 12)
  dir <- tempfile("raefga_"); dir.create(dir); out <- file.path(dir, "r.rae")
  res <- rae(seg, build = "hg38", params = fast(), out = out)
  smp <- read.delim(paste0(out, ".samples.tsv"), check.names = FALSE)
  smp$fga <- NULL; names(smp)[names(smp) == "rae_fga"] <- "cna_burden"
  write.table(smp, paste0(out, ".samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  js <- readLines(paste0(out, ".run.json"))
  writeLines(gsub(', "fga_(threshold|threshold_loss|assembly)": [^,}]+', "", js), paste0(out, ".run.json"))   # as before 3.0.0.9000
  old <- rae_read(out)
  expect_equal(old$samples$rae_fga, res$samples$rae_fga)
  expect_false("cna_burden" %in% names(old$samples))
  expect_false(any(grepl('"fga_', readLines(paste0(out, ".run.json")))))
  expect_true(old$params$fga_assembly)
})

test_that("the command-line tool agrees with rae() for the fga thresholds and --fga-assembly", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  seg <- toy_cohort(n = 12)
  dir <- tempfile("raefgacli_"); dir.create(dir); sf <- file.path(dir, "c.seg")
  write.table(format(seg, scientific = FALSE, trim = TRUE, digits = 15), sf, sep = "\t", quote = FALSE, row.names = FALSE)
  rae(read.delim(sf), build = "hg38", params = fast(fga_threshold = 0.3, fga_threshold_loss = -0.65, fga_assembly = TRUE), out = file.path(dir, "r.rae"))
  system2(cli, c(sf, "--build", "hg38", "--n-perm", "400", "--recleave-every", "100", "--fga-threshold", "0.3", "--fga-threshold-loss", "-0.65",
                 "--fga-assembly", "--out", file.path(dir, "c.rae")), stdout = FALSE, stderr = FALSE)
  expect_identical(readLines(file.path(dir, "c.rae.samples.tsv")), readLines(file.path(dir, "r.rae.samples.tsv")))
})
