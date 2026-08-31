test_that("writePhylip writes a concatenated PHYLIP matrix and a RAxML partition file", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  catdf <- catfullGenes(geneA, geneB, verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".phy")

  writePhylip(catdf, file = tmp, verbose = FALSE)

  lines <- readLines(tmp)
  expect_equal(lines[1], "4 20")
  expect_true(any(grepl("Genus_alpha", lines)))
  expect_true(any(grepl("ACGTACGTACTTGGCCAATT", lines)))

  partition_file <- paste0(sub("[.].*", "", tmp), "_partition_file.txt")
  expect_true(file.exists(partition_file))
  partition_lines <- readLines(partition_file)
  expect_true(any(grepl("DNA, gene1 = 1-10", partition_lines)))
  expect_true(any(grepl("DNA, gene2 = 11-20", partition_lines)))
})

test_that("writePhylip keeps collector identifiers for duplicated accessions (catmultGenes output)", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()
  catdf <- catmultGenes(geneA, geneB,
                        maxspp = TRUE,
                        shortaxlabel = FALSE,
                        missdata = TRUE,
                        verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".phy")

  writePhylip(catdf, file = tmp, verbose = FALSE)

  lines <- readLines(tmp)
  expect_true(any(grepl("Genus_alpha_Coll1", lines)))
  expect_true(any(grepl("Genus_alpha_Coll2", lines)))
  expect_true(any(grepl("^Genus_beta\\s", lines)))
})

test_that("writePhylip genomics = TRUE keeps identifiers for all species, not just duplicates", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()
  catdf <- catmultGenes(geneA, geneB,
                        maxspp = TRUE,
                        shortaxlabel = FALSE,
                        missdata = TRUE,
                        verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".phy")

  writePhylip(catdf, file = tmp, genomics = TRUE, verbose = FALSE)

  lines <- readLines(tmp)
  expect_true(any(grepl("Genus_beta_Coll3", lines)))
  expect_true(any(grepl("Genus_gamma_Coll4", lines)))
})

test_that("writePhylip skips the partition file when partitionfile = FALSE", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  catdf <- catfullGenes(geneA, geneB, verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".phy")

  writePhylip(catdf, file = tmp, partitionfile = FALSE, verbose = FALSE)

  partition_file <- paste0(sub("[.].*", "", tmp), "_partition_file.txt")
  expect_false(file.exists(partition_file))
})

test_that("writePhylip writes only the partition content to `file` when catalignments = FALSE", {
  # When catalignments = FALSE and partitionfile = TRUE, the function writes
  # the partition definitions directly to `file` instead of a separate
  # "*_partition_file.txt" (no concatenated PHYLIP matrix is written at all).
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  catdf <- catfullGenes(geneA, geneB, verbose = FALSE)
  tmp <- withr::local_tempfile(fileext = ".phy")

  writePhylip(catdf, file = tmp, catalignments = FALSE, partitionfile = TRUE, verbose = FALSE)

  expect_true(file.exists(tmp))
  lines <- readLines(tmp)
  expect_true(any(grepl("DNA, gene1 = 1-10", lines)))
})
