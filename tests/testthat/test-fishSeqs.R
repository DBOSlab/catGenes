# Fixture alignment (5 sites) designed so gap-only-column trimming can only
# be verified correctly if it happens AFTER taxon filtering:
#  site 1: A/A/A  - informative, always kept
#  site 2: C/C/C  - informative, always kept
#  site 3: -/-/G  - gap in alpha & beta, but NOT in gamma (which gets fished
#                   out) -> becomes all-gap, and must be trimmed
#  site 4: A/-/A  - gap only in beta -> must be KEPT (not all-gap)
#  site 5: -/-/-  - all-gap regardless -> must be trimmed
.fishseqs_fixture <- function() {
  data.frame(
    species = c("Genus_alpha_V1", "Genus_beta_V2", "Genus_gamma_V3"),
    sequence = c("AC-A-", "AC---", "ACGA-"),
    stringsAsFactors = FALSE
  )
}

test_that("fishSeqs keeps only matching taxa and trims gap-only columns", {
  indir <- withr::local_tempdir()
  # Written directly (not via fastadframe()) to keep exact control over gap
  # placement: fastadframe() would otherwise convert terminal gap runs to "?"
  writeLines(
    c(">Genus_alpha_V1", "AC-A-", ">Genus_beta_V2", "AC---", ">Genus_gamma_V3", "ACGA-"),
    file.path(indir, "gene1.fasta")
  )

  res <- fishSeqs(filepath = indir, taxa = c("alpha", "beta"),
                  save = FALSE, verbose = FALSE)

  expect_equal(names(res), "gene1")
  fished <- res$gene1
  expect_equal(sort(names(fished)), c("Genus_alpha_V1", "Genus_beta_V2"))
  expect_equal(length(fished[[1]]), 3)
  expect_equal(paste(fished[["Genus_alpha_V1"]], collapse = ""), "ACA")
  expect_equal(paste(fished[["Genus_beta_V2"]], collapse = ""), "AC-")
})

test_that("fishSeqs processes every alignment file in the directory", {
  indir <- withr::local_tempdir()
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene1.fasta"))
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene2.fasta"))

  res <- fishSeqs(filepath = indir, taxa = "alpha", save = FALSE, verbose = FALSE)

  expect_equal(sort(names(res)), c("gene1", "gene2"))
})

test_that("fishSeqs saves fished alignments to disk, keeping original format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene1.fasta"))

  fishSeqs(filepath = indir, taxa = "alpha", save = TRUE, dir = outdir, verbose = FALSE)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "gene1.fasta")))
})

test_that("fishSeqs can force a different output format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene1.fasta"))

  fishSeqs(filepath = indir, taxa = "alpha", format = "NEXUS",
          save = TRUE, dir = outdir, verbose = FALSE)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "gene1.nex")))
})

test_that("fishSeqs reads NEXUS and PHYLIP input alignments too", {
  indir <- withr::local_tempdir()
  nexusdframe(.fishseqs_fixture(), file = file.path(indir, "gene1.nex"))
  phylipdframe(.fishseqs_fixture(), file = file.path(indir, "gene2.phy"))

  res <- fishSeqs(filepath = indir, taxa = "alpha", save = FALSE, verbose = FALSE)

  expect_equal(sort(names(res)), c("gene1", "gene2"))
  expect_true("Genus_alpha_V1" %in% names(res$gene1))
  expect_true("Genus_alpha_V1" %in% names(res$gene2))
})

test_that("fishSeqs warns and skips files with none of the requested taxa", {
  indir <- withr::local_tempdir()
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene1.fasta"))

  expect_warning(
    res <- fishSeqs(filepath = indir, taxa = "NoSuchTaxon", save = FALSE, verbose = FALSE),
    "No alignment yielded"
  )
  expect_equal(length(res), 0)
})

test_that("fishSeqs prints progress messages when verbose = TRUE", {
  indir <- withr::local_tempdir()
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene1.fasta"))

  expect_message(
    fishSeqs(filepath = indir, taxa = "alpha", save = FALSE, verbose = TRUE),
    "Kept"
  )
})

test_that("fishSeqs errors when filepath or taxa is missing", {
  expect_error(fishSeqs(filepath = NULL, taxa = "alpha"), "filepath")
  indir <- withr::local_tempdir()
  expect_error(fishSeqs(filepath = indir, taxa = NULL), "taxa")
})

test_that("fishSeqs errors when the directory has no files", {
  indir <- withr::local_tempdir()
  expect_error(fishSeqs(filepath = indir, taxa = "alpha"), "no DNA alignment")
})

test_that("fishSeqs applies rename after fishing and trimming", {
  indir <- withr::local_tempdir()
  fastadframe(.fishseqs_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- c(Genus_alpha_V1 = "Genus_alpha_RENAMED")

  res <- fishSeqs(filepath = indir, taxa = c("alpha", "beta"),
                  rename = lookup, save = FALSE, verbose = FALSE)

  expect_true("Genus_alpha_RENAMED" %in% names(res$gene1))
  expect_true("Genus_beta_V2" %in% names(res$gene1))
  expect_null(attr(res$gene1, "n_renamed"))
})
