test_that("convertAlign converts a FASTA file into NEXUS format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  fastadframe(df, file = file.path(indir, "gene1.fasta"))

  convertAlign(filepath = indir, format = "NEXUS", verbose = FALSE, dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  converted <- list.files(foldername, full.names = TRUE)
  expect_length(converted, 1)
  expect_true(grepl("[.]nex$", converted))
  lines <- readLines(converted)
  expect_true(any(grepl("^#NEXUS", lines)))
})

test_that("convertAlign converts a NEXUS file into PHYLIP format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  nexusdframe(df, file = file.path(indir, "gene1.nex"))

  convertAlign(filepath = indir, format = "PHYLIP", verbose = FALSE, dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  converted <- list.files(foldername, full.names = TRUE)
  expect_length(converted, 1)
  expect_true(grepl("[.]phy$", converted))
})

test_that("convertAlign converts a PHYLIP file into NEXUS format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  phylipdframe(df, file = file.path(indir, "gene1.phy"))

  convertAlign(filepath = indir, format = "NEXUS", verbose = FALSE, dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  converted <- list.files(foldername, full.names = TRUE)
  expect_length(converted, 1)
  lines <- readLines(converted)
  expect_true(any(grepl("Genus_alpha", lines)))
})

test_that("convertAlign converts a NEXUS file into FASTA format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  nexusdframe(df, file = file.path(indir, "gene1.nex"))

  convertAlign(filepath = indir, format = "FASTA", verbose = FALSE, dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  converted <- list.files(foldername, full.names = TRUE)
  expect_length(converted, 1)
  expect_true(grepl("[.]fasta$", converted))
  lines <- readLines(converted)
  expect_equal(lines[1], ">Genus_alpha")
})

test_that("convertAlign skips a file already in the target format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  fastadframe(df, file = file.path(indir, "gene1.fasta"))

  expect_warning(
    expect_message(convertAlign(filepath = indir, format = "FASTA", verbose = FALSE, dir = outdir),
                   "already in FASTA format"),
    "nothing was saved")

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_length(list.files(foldername), 0)
})

test_that("convertAlign rmfiles = TRUE deletes the original input file", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  infile <- file.path(indir, "gene1.fasta")
  fastadframe(df, file = infile)

  convertAlign(filepath = indir, format = "NEXUS", rmfiles = TRUE, verbose = FALSE, dir = outdir)

  expect_false(file.exists(infile))
})

test_that("convertAlign keeps dots in file names so converted files do not collide", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  fastadframe(df, file = file.path(indir, "gene.v1.fasta"))
  fastadframe(df, file = file.path(indir, "gene.v2.fasta"))

  convertAlign(filepath = indir, format = "NEXUS", verbose = FALSE, dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_setequal(basename(list.files(foldername)), c("gene.v1.nex", "gene.v2.nex"))
})

test_that("convertAlign converts FASTA files with wrapped sequences and blank lines", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  writeLines(c(">A", "ACGTACGTAC", "GTACGTACGT", "", ">B", "ACGTACGAAC", "GTACGTACGT", ""),
             file.path(indir, "wrapped.fasta"))

  convertAlign(filepath = indir, format = "PHYLIP", verbose = FALSE, dir = outdir)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  lines <- readLines(file.path(foldername, "wrapped.phy"))
  expect_true(any(grepl("ACGTACGTACGTACGTACGT", lines)))
})

test_that("convertAlign ignores sub-directories and warns when nothing is converted", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  fastadframe(df, file = file.path(indir, "gene1.fasta"))
  dir.create(file.path(indir, "old_results"))

  expect_no_error(
    convertAlign(filepath = indir, format = "NEXUS", verbose = FALSE, dir = outdir)
  )
  expect_warning(
    suppressMessages(
      convertAlign(filepath = indir, format = "FASTA", verbose = FALSE, dir = outdir)
    ),
    "nothing was saved"
  )
})

test_that("convertAlign skips binary and non-alignment files in the input folder", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  fastadframe(df, file = file.path(indir, "gene1.fasta"))
  writeBin(as.raw(c(0x50, 0x4b, 0x03, 0x04, 0x88, 0x39, 0x4c, 0xec, 0x00, 0xff, 0xfd, 0x7c)),
           file.path(indir, "table.xlsx"))
  writeLines(c("species,locus", "A,matK"), file.path(indir, "notes.csv"))

  expect_message(
    convertAlign(filepath = indir, format = "NEXUS", verbose = FALSE, dir = outdir),
    "not a text alignment|not a NEXUS, FASTA or PHYLIP")

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_equal(list.files(foldername), "gene1.nex")
})
