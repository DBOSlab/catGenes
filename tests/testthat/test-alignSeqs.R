# alignSeqs() delegates the actual alignment work to Biostrings::readDNAStringSet()
# and msa::msa()/msaConvert(), which are heavy/slow Bioconductor calls. Both are
# mocked here so tests run fast and don't depend on a working msa installation.

.alignSeqs_mocks <- function() {
  testthat::local_mocked_bindings(readDNAStringSet = function(...) "fake_stringset",
                                  .package = "Biostrings",
                                  .env = parent.frame())
  testthat::local_mocked_bindings(
    msa = function(...) "fake_msa_result",
    msaConvert = function(...) list(nam = c("Genus_alpha", "Genus_beta"),
                                    seq = c("ACGTACGT", "ACGAACGT")),
    .package = "msa",
    .env = parent.frame())
}

test_that("alignSeqs errors when the input directory has no files", {
  empty_dir <- withr::local_tempdir()
  expect_error(alignSeqs(filepath = empty_dir, method = "ClustalW"), "no DNA alignment")
})

test_that("alignSeqs errors when the first file is not FASTA-formatted", {
  indir <- withr::local_tempdir()
  writeLines("not fasta", file.path(indir, "bad.txt"))
  expect_error(alignSeqs(filepath = indir, method = "ClustalW"), 'expected at beginning')
})

test_that("alignSeqs writes a NEXUS-formatted alignment by default", {
  .alignSeqs_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  writeLines(c(">Genus_alpha", "ACGTACGT", ">Genus_beta", "ACGAACGT"),
            file.path(indir, "gene1_ITS.fasta"))

  expect_message(alignSeqs(filepath = indir, method = "ClustalW", dir = outdir),
                 "already saved on disk")

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  nex_files <- list.files(foldername, pattern = "[.]nex$")
  expect_length(nex_files, 1)
})

test_that("alignSeqs writes FASTA output when format = 'FASTA'", {
  .alignSeqs_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  writeLines(c(">Genus_alpha", "ACGTACGT"), file.path(indir, "gene1_ITS.fasta"))

  suppressMessages(alignSeqs(filepath = indir, method = "Muscle", format = "FASTA", dir = outdir))

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_length(list.files(foldername, pattern = "[.]fasta$"), 1)
})

test_that("alignSeqs custom filename prefixes the output file names", {
  .alignSeqs_mocks()
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  writeLines(c(">Genus_alpha", "ACGTACGT"), file.path(indir, "gene1_ITS.fasta"))

  suppressMessages(alignSeqs(filepath = indir, method = "ClustalW",
                             dir = outdir, filename = "myproject"))

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(any(grepl("^myproject_", list.files(foldername))))
})
