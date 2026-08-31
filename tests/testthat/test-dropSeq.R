test_that("dropSeq keeps the accession with the least missing data for a single gene", {
  gene <- fixture_geneA_mult()
  out <- dropSeq(list(A = gene), verbose = FALSE)

  expect_length(out, 1)
  expect_named(out[[1]], c("species", "sequence"))
  # Coll1 has no '?' while Coll2 has 2, so Coll2 must be dropped
  expect_true("Genus_alpha_Coll1" %in% out[[1]]$species)
  expect_false("Genus_alpha_Coll2" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), 3)
})

test_that("dropSeq processes multiple genes with differing numbers of taxa", {
  geneA <- fixture_geneA_mult()
  geneExtra <- fixture_geneA_mult()
  geneExtra[["Extra_taxon_Coll9"]] <- toupper(strsplit("ACGTACGTAC", "")[[1]])

  out <- dropSeq(list(A = geneA, B = geneExtra), verbose = FALSE)

  expect_length(out, 2)
  expect_false("Genus_alpha_Coll2" %in% out[[1]]$species)
  expect_true("Extra_taxon_Coll9" %in% out[[2]]$species)
})

test_that("dropSeq is a no-op when given >= 2 genes with identical taxon counts (documented limitation)", {
  # When more than one gene dataset is supplied and they all have the SAME
  # number of rows, dropSeq's internal `equalnumb()` guard skips the
  # duplicate-removal step entirely - this mirrors the package's own advice to
  # run catmultGenes() first in that scenario.
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- dropSeq(list(A = geneA, B = geneB), verbose = FALSE)

  expect_true("Genus_alpha_Coll2" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), 4)
})

test_that("dropSeq cleanly removes duplicated accessions from catmultGenes() output", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()
  catdf <- catmultGenes(geneA, geneB,
                        maxspp = TRUE,
                        shortaxlabel = FALSE,
                        missdata = TRUE,
                        verbose = FALSE)

  out <- dropSeq(catdf, verbose = FALSE)

  expect_length(out, 2)
  expect_named(out[[1]], c("species", "sequence"))
  expect_false("Genus_alpha_Coll2" %in% out[[1]]$species)
  expect_true("Genus_alpha_Coll1" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})
