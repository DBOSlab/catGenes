# evomodelTest() drives phangorn::modelTest() (slow ML model search) and several
# ape helpers. All of these are mocked so the test runs fast and doesn't depend
# on a working phangorn/ape ML pipeline.

.fake_phy <- function() {
  phy <- list(t1 = 1, t2 = 1, t3 = 1)
  attr(phy, "weight") <- rep(1, 20)
  attr(phy, "type") <- "DNA"
  phy
}

.fake_model_test <- function() {
  data.frame(
    Model = c("JC", "F81", "K80", "HKY", "GTR"),
    AIC = c(120, 115, 110, 100, 105),
    BIC = c(130, 122, 118, 108, 116),
    AICc = c(121, 116, 111, 101, 106),
    logLik = c(-55, -52, -50, -45, -47),
    stringsAsFactors = FALSE
  )
}

.evomodelTest_mocks <- function() {
  testthat::local_mocked_bindings(
    read.phyDat = function(...) .fake_phy(),
    modelTest = function(...) .fake_model_test(),
    pml.control = function(...) list(),
    .package = "phangorn",
    .env = parent.frame())
  testthat::local_mocked_bindings(
    as.DNAbin = function(...) structure(list(), class = "DNAbin"),
    base.freq = function(dna, all = FALSE, ...) {
      if (all) {
        c(a = 0.25, c = 0.25, g = 0.25, t = 0.25, n = 0, `-` = 0)
      } else {
        c(a = 0.25, c = 0.25, g = 0.25, t = 0.25)
      }
    },
    .package = "ape",
    .env = parent.frame())
}

.make_nex_file <- function(dir, name) {
  path <- file.path(dir, name)
  writeLines(c("#NEXUS", "BEGIN DATA;", "MATRIX", "t1 ACGT", "t2 ACGT", "t3 ACGT",
              ";", "END;"), path)
  path
}

test_that("evomodelTest errors when no NEXUS files are found", {
  empty_dir <- withr::local_tempdir()
  expect_error(evomodelTest(nexus.file.path = empty_dir), "No files found")
})

test_that("evomodelTest selects the best model by AIC and writes a report", {
  .evomodelTest_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  .make_nex_file(indir, "geneA.nex")

  res <- evomodelTest(nexus.file.path = indir,
                      model.criteria = "AIC",
                      append.mrbayes.to.nexus = FALSE,
                      verbose = FALSE,
                      dir = outdir)

  expect_equal(res$selected_model, "HKY")
  expect_true(file.exists(res$output_files$report))
  expect_true(file.exists(res$output_files$raw_data))
})

test_that("evomodelTest falls back to AIC with a warning for an invalid criterion", {
  .evomodelTest_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  .make_nex_file(indir, "geneA.nex")

  expect_warning(
    res <- evomodelTest(nexus.file.path = indir,
                        model.criteria = "bogus",
                        append.mrbayes.to.nexus = FALSE,
                        verbose = FALSE,
                        dir = outdir),
    "should be"
  )
  expect_equal(res$selected_model, "HKY")
})

test_that("evomodelTest appends a MrBayes block to a copy of the NEXUS file", {
  .evomodelTest_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  .make_nex_file(indir, "geneA.nex")

  res <- evomodelTest(nexus.file.path = indir,
                      append.mrbayes.to.nexus = TRUE,
                      overwrite.original.nexus = FALSE,
                      verbose = FALSE,
                      dir = outdir)

  expect_true(file.exists(res$output_files$nexus_with_mrbayes))
  lines <- readLines(res$output_files$nexus_with_mrbayes)
  expect_true(any(grepl("begin mrbayes;", lines)))
})

test_that("evomodelTest with multiple partitions writes a combined MrBayes block file", {
  .evomodelTest_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  .make_nex_file(indir, "geneA.nex")
  .make_nex_file(indir, "geneB.nex")

  res <- evomodelTest(nexus.file.path = indir,
                      append.mrbayes.to.nexus = FALSE,
                      verbose = FALSE,
                      dir = outdir)

  expect_named(res, c("individual_files", "combined_analysis"))
  expect_true(file.exists(res$combined_analysis$compact_output_file))
  expect_length(res$combined_analysis$selected_models, 2)
})
