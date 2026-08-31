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

test_that("convertAlign skips a file already in the target format", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  df <- fixture_seq_dframe()
  fastadframe(df, file = file.path(indir, "gene1.fasta"))

  expect_message(convertAlign(filepath = indir, format = "FASTA", verbose = FALSE, dir = outdir),
                 "already in FASTA format")

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
