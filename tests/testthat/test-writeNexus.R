test_that("writeNexus writes an interleaved matrix with a MrBayes block by default", {
  # catfullGenes() derives gene/partition names from the argument expressions
  # (via substitute()), so genes must be passed as short-named variables here.
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  catdf <- catfullGenes(geneA, geneB, verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".nex")

  writeNexus(catdf, file = tmp, verbose = FALSE)

  lines <- readLines(tmp)
  expect_true(any(grepl("^#NEXUS", lines)))
  expect_true(any(grepl("NTAX=4 NCHAR=20", lines)))
  expect_true(any(grepl("INTERLEAVE=YES", lines)))
  expect_true(any(grepl("begin mrbayes;", lines)))
  expect_true(any(grepl("charset geneA = 1-10", lines)))
  expect_true(any(grepl("charset geneB = 11-20", lines)))
})

test_that("writeNexus non-interleaved fully concatenates sequences per taxon", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  catdf <- catfullGenes(geneA, geneB, verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".nex")

  writeNexus(catdf, file = tmp, interleave = FALSE, bayesblock = FALSE, verbose = FALSE)

  lines <- readLines(tmp)
  expect_true(any(grepl("INTERLEAVE=NO", lines)))
  expect_false(any(grepl("begin mrbayes;", lines)))
  expect_true(any(grepl("ACGTACGTACTTGGCCAATT", lines)))
})

test_that("writeNexus errors when given a single non-list gene dataset", {
  tmp <- withr::local_tempfile(fileext = ".nex")
  expect_error(writeNexus(fixture_geneA_full(), file = tmp, verbose = FALSE))
})
