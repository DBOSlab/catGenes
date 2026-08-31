test_that("catfullGenes matches two non-duplicated gene datasets, filling missing data", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()

  out <- catfullGenes(geneA, geneB, verbose = FALSE)

  expect_type(out, "list")
  expect_length(out, 2)
  expect_named(out[[1]], c("species", "sequence"))
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
  expect_setequal(out[[1]]$species, c("Genus_alpha", "Genus_beta", "Genus_gamma", "Genus_delta"))
})

test_that("catfullGenes accepts a single list of gene datasets", {
  genes <- list(A = fixture_geneA_full(), B = fixture_geneB_full())
  out <- catfullGenes(genes, verbose = FALSE)

  expect_length(out, 2)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})

test_that("catfullGenes with missdata = TRUE fills unmatched taxa with '?'", {
  geneA <- fixture_geneA_full()
  geneC <- fixture_geneC_partial()

  out <- catfullGenes(geneA, geneC, missdata = TRUE, verbose = FALSE)

  # geneC lacks Genus_delta, so it should be filled with only '?' characters
  delta_seq <- out[[2]]$sequence[out[[2]]$species == "Genus_delta"]
  expect_true(grepl("^[?]+$", delta_seq))
})

test_that("catfullGenes with missdata = FALSE drops taxa absent from any gene", {
  geneA <- fixture_geneA_full()
  geneC <- fixture_geneC_partial()

  out <- catfullGenes(geneA, geneC, missdata = FALSE, verbose = FALSE)

  # Note: when missdata = FALSE, catfullGenes restores each species' full
  # original identifier (including voucher) regardless of shortaxlabel.
  expect_false(any(grepl("^Genus_delta", out[[1]]$species)))
  expect_setequal(gsub("_V[0-9]$", "", out[[1]]$species),
                  c("Genus_alpha", "Genus_beta", "Genus_gamma"))
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})

test_that("catfullGenes with missdata = FALSE and outgroup keeps the outgroup taxon", {
  geneA <- fixture_geneA_full()
  geneC <- fixture_geneC_partial()

  out <- catfullGenes(geneA, geneC,
                      missdata = FALSE,
                      outgroup = "Genus_delta",
                      verbose = FALSE)

  expect_true(any(grepl("^Genus_delta", out[[1]]$species)))
  # geneC never had delta, so its filled-in sequence should be all '?'
  delta_seq <- out[[2]]$sequence[grepl("^Genus_delta", out[[2]]$species)]
  expect_true(grepl("^[?]+$", delta_seq))
})

test_that("catfullGenes errors when given a single gene dataset (ambiguous with a list-of-genes)", {
  # Passing just one gene alignment cannot be distinguished internally from a
  # single list of >= 2 gene datasets, so it fails downstream rather than with
  # the "at least TWO gene datasets" message.
  expect_error(catfullGenes(fixture_geneA_full(), verbose = FALSE))
})

test_that("catfullGenes errors when a gene dataset has duplicated species", {
  geneA_dup <- fixture_geneA_full()
  geneA_dup[[5]] <- geneA_dup[[1]]
  names(geneA_dup)[5] <- names(geneA_dup)[1]

  expect_error(catfullGenes(geneA_dup, fixture_geneB_full(), verbose = FALSE))
})

test_that("catfullGenes handles a doubtfully-identified (cf.) taxon across genes", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  names(geneA)[2] <- "Genus_cf_beta_V2"
  names(geneB)[2] <- "Genus_cf_beta_V2"

  out <- catfullGenes(geneA, geneB, verbose = FALSE)

  expect_true("Genus_cf_beta" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})

test_that("catfullGenes handles an infraspecific taxon across genes", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()
  names(geneA)[2] <- "Genus_beta_variety_V2"
  names(geneB)[2] <- "Genus_beta_variety_V2"

  out <- catfullGenes(geneA, geneB, verbose = FALSE)

  expect_true("Genus_beta_variety" %in% out[[1]]$species)
  expect_equal(nrow(out[[1]]), nrow(out[[2]]))
})

test_that("catfullGenes shortaxlabel = FALSE keeps voucher identifiers", {
  geneA <- fixture_geneA_full()
  geneB <- fixture_geneB_full()

  out <- catfullGenes(geneA, geneB, shortaxlabel = FALSE, verbose = FALSE)

  expect_true(any(grepl("_V[0-9]$", out[[1]]$species)))
})
