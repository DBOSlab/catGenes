# mineTaxa() performs a GenBank search/download via several rentrez::entrez_*()
# calls and re-parses the result via ape::read.FASTA()/write.FASTA(). All of
# these are mocked below so no network call is ever made.

.mineTaxa_mocks <- function() {
  fake_search <- function(db, term, use_history = FALSE, ...) {
    if (term %in% c("AB111111", "AB222222")) {
      list(ids = list("999"))
    } else {
      list(count = 2, web_history = list(WebEnv = "x", QueryKey = 1))
    }
  }
  fake_fetch <- function(db, id = NULL, web_history = NULL, retmax = NULL,
                         rettype, retmode = NULL, ...) {
    if (rettype == "fasta") {
      paste0(">AB111111.1 Genus species1 voucher1\nACGTACGT\n",
            ">AB222222.1 Genus species2, complete genome\nTTGGCCAA\n")
    } else {
      "LOCUS FAKE\n     /specimen_voucher=\"Cardoso 123\"\n"
    }
  }
  fake_summary <- function(db, id) list(organism = "Genus species1")
  fake_read_fasta <- function(file) {
    out <- list(as.raw(1:4), as.raw(1:4))
    names(out) <- c("AB111111.1 Genus species1 voucher1",
                    "AB222222.1 Genus species2, complete genome")
    class(out) <- "DNAbin"
    out
  }
  fake_write_fasta <- function(...) invisible(NULL)

  testthat::local_mocked_bindings(entrez_search = fake_search,
                                  entrez_fetch = fake_fetch,
                                  entrez_summary = fake_summary,
                                  .package = "rentrez",
                                  .env = parent.frame())
  testthat::local_mocked_bindings(read.FASTA = fake_read_fasta,
                                  write.FASTA = fake_write_fasta,
                                  .package = "ape",
                                  .env = parent.frame())
}

test_that("mineTaxa validates its input arguments before contacting GenBank", {
  expect_error(mineTaxa(term = NULL), "cannot be empty")
  expect_error(mineTaxa(term = ""), "cannot be empty")
  expect_error(mineTaxa(term = "Genus[Organism]", retmax = -1), "retmax must be")
})

test_that("mineTaxa cleans taxon names, adds vouchers, and separates plastomes", {
  .mineTaxa_mocks()
  outdir <- withr::local_tempdir()

  out <- mineTaxa(term = "Genus[Organism]",
                  verbose = FALSE,
                  save = TRUE,
                  dir = outdir,
                  clean.taxa = TRUE,
                  add.voucher = TRUE,
                  plastome.apart = TRUE,
                  rm.duplicated = TRUE)

  expect_true(all(grepl("^Genus_species1_Cardoso123_", names(out))))

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(dir.exists(file.path(foldername, "PLASTOME_GB_FILES")))
  expect_true(any(grepl("AB222222_plastome[.]gb$",
                        list.files(file.path(foldername, "PLASTOME_GB_FILES")))))
})

test_that("mineTaxa without add.voucher omits the voucher segment from names", {
  .mineTaxa_mocks()
  outdir <- withr::local_tempdir()

  out <- mineTaxa(term = "Genus[Organism]",
                  verbose = FALSE,
                  save = TRUE,
                  dir = outdir,
                  clean.taxa = TRUE,
                  add.voucher = FALSE,
                  plastome.apart = FALSE)

  expect_true(all(grepl("^Genus_species1_AB", names(out))))
})

test_that("mineTaxa clean.taxa = FALSE returns the raw FASTA text", {
  .mineTaxa_mocks()
  outdir <- withr::local_tempdir()

  out <- mineTaxa(term = "Genus[Organism]",
                  verbose = FALSE,
                  save = TRUE,
                  dir = outdir,
                  clean.taxa = FALSE)

  expect_true(is.character(out))
  expect_true(grepl("^>AB111111", out))
})

test_that(".rm_duplicated_seqs keeps only the longest sequence per species", {
  # Species-name extraction inside .rm_duplicated_seqs assumes a Latin-binomial
  # shape (letters only, no digits) for the specific epithet.
  seqs <- list(short = as.raw(1:2), long = as.raw(1:5), other = as.raw(1:3))
  names(seqs) <- c("Genus_alpha_v1", "Genus_alpha_v2", "Genus_beta_v1")

  out <- catGenes:::.rm_duplicated_seqs(seqs)

  expect_equal(names(out), c("Genus_alpha_v2", "Genus_beta_v1"))
})

test_that(".rm_duplicated_seqs is a no-op on an empty list", {
  expect_equal(catGenes:::.rm_duplicated_seqs(list()), list())
})
