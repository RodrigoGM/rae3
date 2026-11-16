# Plotting functions and the germline switch. toy_cohort() and fast() are in helper-toy.R.

res_with_normals <- function(build = "custom", ...) {
  seg <- toy_cohort()
  nrm <- seg[seg$ID %in% sprintf("S%02d", 1:6), ]; nrm$ID <- sub("S", "N", nrm$ID)
  rae(seg, normals = nrm, build = build, params = fast(...), cnv = "none")
}

test_that("germline overlaps are flagged by default and moved out with filter_germline_overlaps = TRUE", {
  a <- res_with_normals()
  expect_true("germline_overlap" %in% names(a$roi))
  g <- a$roi$germline_overlap
  expect_true(any(g))                                  # the chr2 loss is present in the normals
  expect_null(a$roi_germline)
  seg <- toy_cohort(); nrm <- seg[seg$ID %in% sprintf("S%02d", 1:6), ]; nrm$ID <- sub("S", "N", nrm$ID)
  b <- rae(seg, normals = nrm, build = "custom", params = fast(), cnv = "none", filter_germline_overlaps = TRUE)
  expect_false(any(b$roi$germline_overlap))
  expect_equal(nrow(b$roi) + nrow(b$roi_germline), nrow(a$roi))
  expect_identical(b$regions, a$regions)               # q-values are never changed
  # switching back restores the full table
  expect_equal(rae_germline(b, filter = FALSE)$roi[, names(a$roi)], a$roi, ignore_attr = TRUE)
})

test_that("the density plot is recomputed on the run's own stream (E reproduced)", {
  r <- rae(toy_cohort(), build = "custom", params = fast())
  for (i in c(1, 7, 40)) {
    s <- r$samples[i, ]; sg <- r$segments[r$segments$ID == s$ID, ]
    f <- rae_parameterize(sg$seg.mean.centred + s$shift, sg$seg.sd, sg$num.mark, dn = s$dn, params = r$params, stream = i - 1)
    expect_equal(unname(f$E[c(1, 3)]), c(s$E_A0, s$E_D0), tolerance = 1e-5)
  }
})

test_that("every plotting function draws on png() and cairo_pdf()", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("ggrepel")
  r <- res_with_normals("hg38", arm_test = TRUE, arm_n_perm = 200L)
  expect_equal(r$build, "hg38")
  expect_false(is.null(r$arms))
  expect_true(nrow(r$diagnostics$hotspots) > 0)
  genes <- data.frame(gene = c("GENEA", "GENEB"), chrom = c("2", "3"), start = c(16e6, 20.5e6), end = c(16.1e6, 20.6e6))
  plots <- list(
    rae_plot_genome(r, genes = genes), rae_plot_genome(r, chrom = c(2, "3q"), genes = NULL, show = "all"),
    rae_plot_genome(rae_germline(r, TRUE), genes = genes),
    rae_plot_locus(r, 2, 10e6, 25e6, genes = genes), rae_plot_locus(r, 3, direction = "gain", genes = NULL),
    rae_plot_noise(r), rae_plot_noise(r, against = "E"), rae_plot_density(r, 1), rae_plot_sample(r, "S02"),
    rae_plot_E(r), rae_plot_diploid_offset(r), rae_plot_segments(r, 2, 10e6, 25e6),
    rae_plot_hotspots(r, "1p"), rae_plot_region_sizes(r), rae_plot_arms(r), rae_plot_arms(r, "bar"))
  for (p in plots) expect_s3_class(p, "ggplot")
  r0 <- r; r0$arms <- NULL
  expect_error(rae_plot_arms(r0), "arm_test")
  d <- withr::local_tempdir()
  for (k in seq_along(plots)) for (e in c("png", "pdf")) {
    f <- file.path(d, paste0(k, ".", e)); rae_save(plots[[k]], f, 100, 70)
    expect_gt(file.size(f), 1000)
  }
  skip_if_not_installed("ComplexHeatmap"); skip_if_not_installed("vegan")
  h <- rae_heatmap(r); expect_s4_class(h, "Heatmap")
  f <- file.path(d, "h.png"); rae_save(h, f, 100, 80); expect_gt(file.size(f), 1000)
})

test_that("genome labels: focal genes, otherwise cytobands, germline mark", {
  cb <- rae_cytobands("hg38")
  ev <- data.frame(chrom = c(10, 10, 8), start = c(87.8e6, 80e6, 1e6), end = c(88.0e6, 95e6, 2e6),
                   size_bp = c(2e5, 15e6, 1e6), germline_overlap = c(FALSE, FALSE, TRUE))
  gt <- data.frame(gene = c("PTEN", "X1"), chrom = c("10", "10"), start = c(87.86e6, 90e6), end = c(87.97e6, 90.1e6))
  lab <- rae3:::.event_label(ev, cb, gt, 3e6, 10, "Г")
  expect_equal(lab[1], "PTEN")
  expect_match(lab[2], "^10q2")                       # broad: cytoband span
  expect_equal(lab[3], "Г")
})

test_that("cytobands are bundled for the three builds", {
  for (b in c("hg19", "hg38", "chm13")) {
    cb <- rae_cytobands(b)
    expect_equal(nrow(cb), 862); expect_true(all(c("acen") %in% cb$stain))
  }
  expect_null(rae_cytobands("custom"))
})
