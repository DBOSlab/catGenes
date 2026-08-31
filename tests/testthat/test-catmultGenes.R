test_that("catmultGenes matches multi-accession genes, keeping all accessions of duplicated species", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB, verbose = FALSE)

  expect_length(out, 2)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
  expect_equal(sum(grepl("^Genus_alpha", out[[1]]$species)), 2)
})

test_that("catmultGenes maxspp = TRUE collapses non-duplicated species to short labels", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB, maxspp = TRUE, verbose = FALSE)

  expect_true("Genus_beta" %in% out[[1]]$species)
  expect_true("Genus_gamma" %in% out[[1]]$species)
})

test_that("catmultGenes maxspp = FALSE keeps collector identifiers for all species", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB, maxspp = FALSE, verbose = FALSE)

  expect_true("Genus_beta_Coll3" %in% out[[1]]$species)
  expect_true("Genus_gamma_Coll4" %in% out[[1]]$species)
})

test_that("catmultGenes errors when no species is duplicated with multiple accessions", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()

  expect_error(catmultGenes(geneA, geneB, verbose = FALSE), "do not include species duplicated")
})

test_that("catmultGenes handles a doubtfully-identified (aff.) accession", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()
  names(geneA)[3] <- "Genus_aff_beta_Coll3"
  names(geneB)[3] <- "Genus_aff_beta_Coll3"

  out <- catmultGenes(geneA, geneB, verbose = FALSE)

  expect_true("Genus_aff_beta" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})

test_that("catmultGenes shortaxlabel = FALSE keeps collector identifiers throughout", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB, shortaxlabel = FALSE, verbose = FALSE)

  expect_true(all(grepl("_Coll[0-9]$", out[[1]]$species)))
})
