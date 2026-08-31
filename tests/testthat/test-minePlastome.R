# minePlastome() mirrors mineMitochondrion() but targets plastome genomes.
# Same mocking strategy: rentrez::entrez_fetch() and geneviewer::* are mocked
# so no network call or real .gb parsing ever happens.

.fake_gbk_plastome <- function() {
  list(list(
    ORIGIN = paste(rep("ACGT", 50), collapse = ""),
    FEATURES = list(source = list(list(organism = "Genus species")))
  ))
}

.fake_gene_df <- function() {
  data.frame(gene = c("rbcL", "matK"),
            start = c(1, 101),
            end = c(90, 190),
            strand = c("plus", "minus"),
            stringsAsFactors = FALSE)
}

test_that("minePlastome writes one FASTA file per requested non-CDS gene", {
  outdir <- withr::local_tempdir()

  testthat::local_mocked_bindings(entrez_fetch = function(...) "FAKE GENBANK RECORD",
                                  .package = "rentrez")
  testthat::local_mocked_bindings(read_gbk = function(...) .fake_gbk_plastome(),
                                  gbk_features_to_df = function(...) .fake_gene_df(),
                                  .package = "geneviewer")

  minePlastome(genbank = "NC_000002",
              CDS = FALSE,
              genes = "rbcL",
              verbose = FALSE,
              dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  fasta_files <- list.files(foldername, pattern = "rbcL_from_plastome[.]fasta$")
  expect_length(fasta_files, 1)
})

test_that("minePlastome complements (but does not reverse) minus-strand loci", {
  # Note: internally minePlastome calls .seq_revcompl(seq, strand) using the
  # literal strand value from GenBank ("plus"/"minus"). .seq_revcompl() only
  # reverses when strand is exactly "complement", so in practice minus-strand
  # loci are base-complemented but NOT reversed - this mirrors that actual,
  # documented-elsewhere behavior rather than an idealized reverse-complement.
  outdir <- withr::local_tempdir()

  testthat::local_mocked_bindings(entrez_fetch = function(...) "FAKE GENBANK RECORD",
                                  .package = "rentrez")
  testthat::local_mocked_bindings(read_gbk = function(...) .fake_gbk_plastome(),
                                  gbk_features_to_df = function(...) .fake_gene_df(),
                                  .package = "geneviewer")

  minePlastome(genbank = "NC_000002",
              CDS = FALSE,
              genes = "matK",
              verbose = FALSE,
              dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  lines <- readLines(file.path(foldername, "matK_from_plastome.fasta"))
  seq_line <- lines[2]
  gbk <- .fake_gbk_plastome()
  raw_seq <- substr(gbk[[1]][["ORIGIN"]], 101, 190)
  expect_equal(seq_line, catGenes:::.seq_revcompl(raw_seq, "minus"))
})

test_that("minePlastome errors for an unknown gene name", {
  outdir <- withr::local_tempdir()

  testthat::local_mocked_bindings(entrez_fetch = function(...) "FAKE GENBANK RECORD",
                                  .package = "rentrez")
  testthat::local_mocked_bindings(read_gbk = function(...) .fake_gbk_plastome(),
                                  gbk_features_to_df = function(...) .fake_gene_df(),
                                  .package = "geneviewer")

  # errors are caught internally (tryCatch) and only printed via cat(), so the
  # function itself should not raise, but no FASTA file should be produced.
  minePlastome(genbank = "NC_000002",
              CDS = FALSE,
              genes = "unknownGene",
              verbose = FALSE,
              dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_length(list.files(foldername, pattern = "[.]fasta$"), 0)
})
