.renametaxa_fixture <- function() {
  data.frame(
    species = c("FC520", "FC459", "Genus_gamma_V3"),
    sequence = c("ACGTACGTAC", "ACGAACGTAC", "ACGTTCGTAC"),
    stringsAsFactors = FALSE
  )
}

test_that("renameTaxa renames matched sequences and leaves the rest unchanged", {
  indir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- c(FC520 = "Moldenhawera_blanchetiana_FC520",
             FC459 = "Tachigali_costaricensis_FC459")

  res <- renameTaxa(filepath = indir, lookup = lookup, save = FALSE, verbose = FALSE)

  expect_equal(
    sort(names(res$gene1)),
    sort(c("Moldenhawera_blanchetiana_FC520", "Tachigali_costaricensis_FC459",
          "Genus_gamma_V3"))
  )
})

test_that("renameTaxa accepts a two-column data.frame lookup", {
  indir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- data.frame(
    old = c("FC520", "FC459"),
    new = c("Moldenhawera_blanchetiana_FC520", "Tachigali_costaricensis_FC459"),
    stringsAsFactors = FALSE
  )

  res <- renameTaxa(filepath = indir, lookup = lookup, save = FALSE, verbose = FALSE)

  expect_true("Moldenhawera_blanchetiana_FC520" %in% names(res$gene1))
})

test_that("renameTaxa saves renamed alignments to a new folder by default", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- c(FC520 = "Moldenhawera_blanchetiana_FC520")

  renameTaxa(filepath = indir, lookup = lookup, dir = outdir, verbose = FALSE)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "gene1.fasta")))

  # original file must remain untouched
  original <- ape::read.FASTA(file.path(indir, "gene1.fasta"))
  expect_true("FC520" %in% names(original))
})

test_that("renameTaxa can overwrite the original files in place", {
  indir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- c(FC520 = "Moldenhawera_blanchetiana_FC520")

  renameTaxa(filepath = indir, lookup = lookup, overwrite = TRUE, verbose = FALSE)

  overwritten <- ape::read.FASTA(file.path(indir, "gene1.fasta"))
  expect_true("Moldenhawera_blanchetiana_FC520" %in% names(overwritten))
  expect_false("FC520" %in% names(overwritten))
})

test_that("renameTaxa can force a different output format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- c(FC520 = "Moldenhawera_blanchetiana_FC520")

  renameTaxa(filepath = indir, lookup = lookup, format = "NEXUS",
            dir = outdir, verbose = FALSE)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "gene1.nex")))
})

test_that("renameTaxa prints progress messages when verbose = TRUE", {
  indir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  lookup <- c(FC520 = "Moldenhawera_blanchetiana_FC520")

  expect_message(
    renameTaxa(filepath = indir, lookup = lookup, save = FALSE, verbose = TRUE),
    "Renamed"
  )
})

test_that("renameTaxa errors when filepath or lookup is missing", {
  expect_error(renameTaxa(filepath = NULL, lookup = c(a = "b")), "filepath")
  indir <- withr::local_tempdir()
  expect_error(renameTaxa(filepath = indir, lookup = NULL), "lookup")
})

test_that("renameTaxa errors when the directory has no files", {
  indir <- withr::local_tempdir()
  expect_error(renameTaxa(filepath = indir, lookup = c(a = "b")), "no DNA alignment")
})

test_that("renameTaxa errors when lookup has no names", {
  indir <- withr::local_tempdir()
  fastadframe(.renametaxa_fixture(), file = file.path(indir, "gene1.fasta"))

  expect_error(
    renameTaxa(filepath = indir, lookup = c("a", "b"), save = FALSE, verbose = FALSE),
    "named vector"
  )
})
