# rae_read(): the command-line tool's output read back into the object rae() returns.

tables <- c("regions", "roi", "samples", "segments", "calls")

test_that("the output set has eight files, including the null-model internals", {
  res <- rae(toy_cohort(), build = "hg19", params = fast())
  expect_length(res$files, 8)
  expect_true(any(grepl("\\.null\\.tsv$", res$files)))
  nl <- utils::read.delim(grep("\\.null\\.tsv$", res$files, value = TRUE))
  expect_named(nl, c("item", "chrom", "value"))
  expect_setequal(unique(nl$item), c("chrom_len", "cen_start", "cen_end", "hotspot", "sampled_cut", "ubp_bp", "cleaved_bp"))
})

test_that("run.json is valid JSON: no bare NA, full-precision parameters", {
  res <- rae(toy_cohort()[, 1:6], build = "custom", params = fast())   # virtual grid: grid_spacing is set
  js <- paste(res$manifest, collapse = "\n")
  expect_false(grepl(":\\s*NA\\b", js))
  expect_match(js, "\"beta_a1\": 1.4426950408889634")
  skip_if_not_installed("jsonlite")
  expect_silent(jsonlite::fromJSON(js))
})

test_that("rae_read() rebuilds the object rae() returned", {
  res <- rae(toy_cohort(), build = "hg19", params = fast(), out = tempfile("rd_"))
  prefix <- sub("\\.regions\\.tsv$", "", grep("regions", res$files, value = TRUE))
  rd <- rae_read(prefix)
  expect_s3_class(rd, "rae_result")
  for (t in tables) expect_identical(rd[[t]], res[[t]], info = t)
  expect_identical(rd$manifest, res$manifest)
  expect_identical(rd$build, res$build)
  expect_identical(rd$files, res$files)
  expect_equal(rd$params, res$params)
  expect_equal(rd$diagnostics, res$diagnostics)
  expect_true(rd$diagnostics$has_hotspots)
  expect_gt(nrow(rd$diagnostics$hotspots), 0)
})

test_that("rae_read() applies the germline switch like rae()", {
  seg <- toy_cohort()
  f <- tempfile(fileext = ".tsv")
  writeLines(c("chrom\tstart\tend\ttype\taf", "2\t15000001\t16000000\tloss\t0.05", "2\t15500001\t17000000\tloss\t0.02"), f)
  res <- rae(seg, build = "custom", params = fast(), cnv = f, filter_germline_overlaps = TRUE)
  prefix <- sub("\\.regions\\.tsv$", "", grep("regions", res$files, value = TRUE))
  rd <- rae_read(prefix, filter_germline_overlaps = TRUE)
  expect_identical(rd$roi, res$roi)
  expect_identical(rd$roi_germline, res$roi_germline)
  expect_gt(nrow(rd$roi_germline), 0)
})

test_that("rae_read() accepts an output without the null-model file and says what is missing", {
  res <- rae(toy_cohort(), build = "custom", params = fast(), out = tempfile("old_"))
  prefix <- sub("\\.regions\\.tsv$", "", grep("regions", res$files, value = TRUE))
  file.remove(paste0(prefix, ".null.tsv"))
  expect_message(rd <- rae_read(prefix), "null.tsv")
  expect_identical(rd$regions, res$regions)
  expect_null(rd$diagnostics$ubp_bp)
  expect_error(rae_plot_region_sizes(rd), "null-model")
  expect_error(rae_read(tempfile("missing_")), "not found")
})

test_that("rae_read() on the command-line tool's output equals the R result", {
  cli <- Sys.getenv("RAE_CLI", "")
  skip_if(cli == "" || !file.exists(cli), "RAE_CLI not set")
  seg <- toy_cohort()
  d <- tempfile("raecli_"); dir.create(d)
  f <- file.path(d, "toy.seg")
  write.table(format(seg, scientific = FALSE, trim = TRUE, digits = 15), f, sep = "\t", quote = FALSE, row.names = FALSE)
  res <- rae(read.delim(f), build = "hg19", params = fast(), out = file.path(d, "r.rae"))
  system2(cli, c(f, "--build", "hg19", "--n-perm", "400", "--recleave-every", "100", "--out", file.path(d, "c.rae")),
          stdout = FALSE, stderr = FALSE)
  rd <- rae_read(file.path(d, "c.rae"))
  for (t in tables) expect_identical(rd[[t]], res[[t]], info = t)
  expect_equal(rd$params, res$params)
  expect_equal(rd$diagnostics, res$diagnostics)
  skip_if_not_installed("ggplot2"); skip_if_not_installed("ggrepel")
  expect_s3_class(rae_plot_genome(rd, genes = NULL), "ggplot")
  expect_s3_class(rae_plot_region_sizes(rd), "ggplot")
  expect_s3_class(rae_plot_hotspots(rd, "1p"), "ggplot")
})
