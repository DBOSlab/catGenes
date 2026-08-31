test_that("nexusdframe writes a well-formed NEXUS file from a data.frame", {
  df <- fixture_seq_dframe()
  tmp <- withr::local_tempfile(fileext = ".nex")

  nexusdframe(df, file = tmp)

  expect_true(file.exists(tmp))
  lines <- readLines(tmp)
  expect_true(any(grepl("^#NEXUS", lines)))
  expect_true(any(grepl("BEGIN DATA;", lines)))
  expect_true(any(grepl("MATRIX", lines)))
  expect_true(any(grepl("Genus_alpha", lines)))
  expect_true(any(grepl(paste0("NTAX=", nrow(df)), lines)))
})

test_that("nexusdframe accepts a list-formatted NEXUS object (as from ape::read.nexus.data)", {
  gene <- fixture_geneA_full()
  tmp <- withr::local_tempfile(fileext = ".nex")

  nexusdframe(gene, file = tmp)

  lines <- readLines(tmp)
  expect_true(any(grepl(paste0("NTAX=", length(gene)), lines)))
  expect_true(any(grepl("Genus_alpha_V1", lines)))
})

test_that("nexusdframe errors when given a single-column data.frame", {
  df <- data.frame(species = c("a", "b"))
  tmp <- withr::local_tempfile(fileext = ".nex")
  expect_error(nexusdframe(df, file = tmp), "two-column")
})

test_that("nexusdframe dropmisseq removes fully-missing sequences", {
  df <- data.frame(species = c("sp1", "sp2"),
                   sequence = c("????????", "ACGTACGT"),
                   stringsAsFactors = FALSE)
  tmp <- withr::local_tempfile(fileext = ".nex")

  nexusdframe(df, file = tmp, dropmisseq = TRUE)

  lines <- readLines(tmp)
  expect_false(any(grepl("sp1", lines)))
  expect_true(any(grepl("sp2", lines)))
})

test_that("phylipdframe writes a well-formed PHYLIP file", {
  df <- fixture_seq_dframe()
  tmp <- withr::local_tempfile(fileext = ".phy")

  phylipdframe(df, file = tmp)

  lines <- readLines(tmp)
  expect_true(grepl(paste0("^", nrow(df), " "), lines[1]))
  expect_true(any(grepl("Genus_alpha", lines)))
})

test_that("phylipdframe errors when given a single-column data.frame", {
  df <- data.frame(species = c("a", "b"))
  tmp <- withr::local_tempfile(fileext = ".phy")
  expect_error(phylipdframe(df, file = tmp), "two-column")
})

test_that("fastadframe writes a well-formed FASTA file", {
  df <- fixture_seq_dframe()
  tmp <- withr::local_tempfile(fileext = ".fasta")

  fastadframe(df, file = tmp)

  lines <- readLines(tmp)
  expect_equal(lines[1], ">Genus_alpha")
  expect_equal(lines[2], "ACGTACGTAC")
  expect_equal(sum(grepl("^>", lines)), nrow(df))
})

test_that("fastadframe drops sequences that are entirely missing data", {
  df <- data.frame(species = c("sp1", "sp2"),
                   sequence = c("????", "ACGT"),
                   stringsAsFactors = FALSE)
  tmp <- withr::local_tempfile(fileext = ".fasta")

  fastadframe(df, file = tmp, dropmisseq = TRUE)

  lines <- readLines(tmp)
  expect_false(any(grepl("sp1", lines)))
  expect_true(any(grepl("sp2", lines)))
})

test_that("fastadframe keeps fully-missing sequences when dropmisseq = FALSE", {
  df <- data.frame(species = c("sp1", "sp2"),
                   sequence = c("????", "ACGT"),
                   stringsAsFactors = FALSE)
  tmp <- withr::local_tempfile(fileext = ".fasta")

  fastadframe(df, file = tmp, dropmisseq = FALSE)

  lines <- readLines(tmp)
  expect_true(any(grepl("sp1", lines)))
})
