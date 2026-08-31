# mineMitochondrion() connects to GenBank via rentrez::entrez_fetch() and parses
# the record via geneviewer::read_gbk()/gbk_features_to_df(). All three are
# mocked below so no network call or real .gb parsing ever happens.

.fake_gbk <- function() {
  list(list(
    ORIGIN = paste(rep("ACGT", 50), collapse = ""),
    FEATURES = list(source = list(list(organism = "Genus species")))
  ))
}

.fake_cds_df <- function() {
  data.frame(gene = c("COX1", "ND4L"),
            start = c(1, 101),
            end = c(90, 190),
            strand = c("plus", "plus"),
            stringsAsFactors = FALSE)
}

test_that("mineMitochondrion writes one FASTA file per requested gene", {
  outdir <- withr::local_tempdir()

  testthat::local_mocked_bindings(entrez_fetch = function(...) "FAKE GENBANK RECORD",
                                  .package = "rentrez")
  testthat::local_mocked_bindings(read_gbk = function(...) .fake_gbk(),
                                  gbk_features_to_df = function(...) .fake_cds_df(),
                                  .package = "geneviewer")

  mineMitochondrion(genbank = "NC_000001",
                    CDS = TRUE,
                    genes = "COX1",
                    verbose = FALSE,
                    dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  fasta_files <- list.files(foldername, pattern = "COX1_from_mitochondrion[.]fasta$")
  expect_length(fasta_files, 1)

  lines <- readLines(file.path(foldername, fasta_files))
  expect_true(any(grepl("^>", lines)))
})

test_that("mineMitochondrion uses the supplied taxon/voucher in the FASTA header", {
  outdir <- withr::local_tempdir()

  testthat::local_mocked_bindings(entrez_fetch = function(...) "FAKE GENBANK RECORD",
                                  .package = "rentrez")
  testthat::local_mocked_bindings(read_gbk = function(...) .fake_gbk(),
                                  gbk_features_to_df = function(...) .fake_cds_df(),
                                  .package = "geneviewer")

  mineMitochondrion(genbank = "NC_000001",
                    taxon = "My_taxon",
                    voucher = "Coll123",
                    CDS = TRUE,
                    genes = "COX1",
                    verbose = FALSE,
                    dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  lines <- readLines(file.path(foldername, "COX1_from_mitochondrion.fasta"))
  expect_true(any(grepl("^>My_taxon_Coll123_", lines)))
})

test_that("mineMitochondrion keeps the .gb file by default and removes it when rm_gb_files = TRUE", {
  outdir <- withr::local_tempdir()

  testthat::local_mocked_bindings(entrez_fetch = function(...) "FAKE GENBANK RECORD",
                                  .package = "rentrez")
  testthat::local_mocked_bindings(read_gbk = function(...) .fake_gbk(),
                                  gbk_features_to_df = function(...) .fake_cds_df(),
                                  .package = "geneviewer")

  mineMitochondrion(genbank = "NC_000001",
                    CDS = TRUE,
                    genes = "COX1",
                    rm_gb_files = TRUE,
                    verbose = FALSE,
                    dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_length(list.files(foldername, pattern = "[.]gb$"), 0)
})
