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

test_that("catmultGenes missdata = FALSE fully matches genes without an outgroup", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB, missdata = FALSE, verbose = FALSE)

  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
  expect_equal(sum(grepl("^Genus_alpha", out[[1]]$species)), 2)
})

test_that("catmultGenes missdata = FALSE with an outgroup keeps it separate from the ingroup match", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB,
                      missdata = FALSE,
                      outgroup = "Genus_gamma",
                      verbose = FALSE)

  expect_true("Genus_gamma" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})

test_that("catmultGenes errors when infraspecific ID is applied inconsistently across accessions", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()
  names(geneA)[3] <- "Genus_beta_variety_Coll3"
  names(geneB)[3] <- "Genus_beta_Coll5"

  expect_error(catmultGenes(geneA, geneB, verbose = FALSE),
              "identified at infraspecific level")
})

test_that("catmultGenes maxspp = TRUE with shortaxlabel = FALSE restores original identifiers", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()

  out <- catmultGenes(geneA, geneB, maxspp = TRUE, shortaxlabel = FALSE, verbose = FALSE)

  expect_true(all(grepl("_Coll[0-9]$", out[[1]]$species)))
})

test_that("catmultGenes verbose = TRUE reports progress across 3+ gene datasets", {
  geneA <- fixture_geneA_mult()
  geneB <- fixture_geneB_mult()
  geneC <- fixture_geneA_mult()
  names(geneC) <- paste0(names(geneC), "b")

  expect_message(catmultGenes(geneA, geneB, geneC, verbose = TRUE),
                 "Full gene match is finished")
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
