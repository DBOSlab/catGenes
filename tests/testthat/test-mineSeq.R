# mineSeq() downloads from GenBank via ape::read.GenBank(); mock it so no
# network call is ever made.

.fake_read_genbank <- function(access.nb, species.names = TRUE, as.character = FALSE) {
  out <- as.list(paste0("fake_seq_for_", access.nb))
  names(out) <- access.nb
  attr(out, "species") <- paste0("Genus_species", seq_along(access.nb))
  attr(out, "description") <- paste0("desc_", access.nb)
  out
}

test_that("mineSeq downloads and names sequences using Species/Voucher columns", {
  indir <- withr::local_tempdir()
  df <- data.frame(Species = c("Genus species1", "Genus species2"),
                   Voucher = c("Cardoso 1", "Cardoso 2"),
                   ITS = c("MN111111", "MN222222"),
                   stringsAsFactors = FALSE)

  testthat::local_mocked_bindings(read.GenBank = .fake_read_genbank, .package = "ape")

  out <- mineSeq(inputdf = df, gb.colnames = "ITS", verbose = FALSE, save = FALSE)

  expect_named(out, "ITS")
  expect_true(any(grepl("^Genus_species1_Cardoso1_MN111111$", names(out$ITS))))
})

test_that("mineSeq works without Species/Voucher columns, using GenBank taxonomy", {
  df <- data.frame(ITS = c("MN111111", "MN222222"), stringsAsFactors = FALSE)

  testthat::local_mocked_bindings(read.GenBank = .fake_read_genbank, .package = "ape")

  out <- mineSeq(inputdf = df, gb.colnames = "ITS", verbose = FALSE, save = FALSE)

  expect_true(any(grepl("^Genus_species1_MN111111$", names(out$ITS))))
})

test_that("mineSeq with Species but no Voucher column omits the voucher segment", {
  df <- data.frame(Species = c("Genus species1", "Genus species2"),
                   ITS = c("MN111111", "MN222222"),
                   stringsAsFactors = FALSE)

  testthat::local_mocked_bindings(read.GenBank = .fake_read_genbank, .package = "ape")

  out <- mineSeq(inputdf = df, gb.colnames = "ITS", verbose = FALSE, save = FALSE)

  expect_true(any(grepl("^Genus_species1_MN111111$", names(out$ITS))))
})

test_that("mineSeq with Voucher but no Species column falls back to GenBank taxonomy", {
  df <- data.frame(Voucher = c("Cardoso 1", "Cardoso 2"),
                   ITS = c("MN111111", "MN222222"),
                   stringsAsFactors = FALSE)

  testthat::local_mocked_bindings(read.GenBank = .fake_read_genbank, .package = "ape")

  out <- mineSeq(inputdf = df, gb.colnames = "ITS", verbose = FALSE, save = FALSE)

  expect_true(any(grepl("^Genus_species1_Cardoso1_MN111111$", names(out$ITS))))
})

test_that("mineSeq as.character = TRUE is passed through to ape::read.GenBank", {
  captured_as_character <- NULL
  fake_read <- function(access.nb, species.names = TRUE, as.character = FALSE) {
    captured_as_character <<- as.character
    out <- as.list(paste0("fake_seq_for_", access.nb))
    names(out) <- access.nb
    attr(out, "species") <- paste0("Genus_species", seq_along(access.nb))
    out
  }
  df <- data.frame(ITS = "MN111111", stringsAsFactors = FALSE)

  testthat::local_mocked_bindings(read.GenBank = fake_read, .package = "ape")

  mineSeq(inputdf = df, gb.colnames = "ITS", as.character = TRUE, verbose = FALSE, save = FALSE)

  expect_true(captured_as_character)
})

test_that("mineSeq saves a FASTA file per gene when save = TRUE", {
  outdir <- withr::local_tempdir()
  df <- data.frame(Species = "Genus species1",
                   Voucher = "Cardoso 1",
                   ITS = "MN111111",
                   stringsAsFactors = FALSE)

  testthat::local_mocked_bindings(read.GenBank = .fake_read_genbank, .package = "ape")
  testthat::local_mocked_bindings(write.dna = function(...) invisible(NULL), .package = "ape")

  out <- mineSeq(inputdf = df, gb.colnames = "ITS", verbose = FALSE,
                save = TRUE, dir = outdir, filename = "test_seqs")

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(dir.exists(foldername))
})

test_that("mineSeq errors when GenBank accessions are duplicated", {
  df <- data.frame(ITS = c("MN111111", "MN111111"), stringsAsFactors = FALSE)
  expect_error(mineSeq(inputdf = df, gb.colnames = "ITS", verbose = FALSE, save = FALSE),
              "duplicated")
})
