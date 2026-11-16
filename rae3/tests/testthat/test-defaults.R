# Default output must not change when opt-in features are added. The fixtures were produced by
# tools/make_default_fixtures.R from commit e62b528.

skip_if_not_installed("jsonlite")
fixture_dir <- function(nm) testthat::test_path("fixtures", "defaults", nm)

test_that("default outputs are byte-identical to the e62b528 fixtures", {
  for (nm in names(default_cases())) {
    dir <- run_case(default_cases()[[nm]])
    for (s in c("regions", "roi", "samples", "segments", "calls", "null"))
      expect_identical(readLines(file.path(dir, paste0("toy.", s, ".tsv"))),
                       readLines(file.path(fixture_dir(nm), paste0(s, ".tsv"))),
                       label = paste(nm, s))
  }
})

test_that("run.json differs only by new parameter keys", {
  for (nm in names(default_cases())) {
    dir <- run_case(default_cases()[[nm]])
    new <- jsonlite::fromJSON(file.path(dir, "toy.run.json"), simplifyVector = FALSE)
    old <- jsonlite::fromJSON(file.path(fixture_dir(nm), "run.json"), simplifyVector = FALSE)
    new$parameters <- new$parameters[names(new$parameters) %in% names(old$parameters)]
    drop <- c("input", "version")
    expect_identical(new[setdiff(names(new), drop)], old[setdiff(names(old), drop)], label = nm)
  }
})
