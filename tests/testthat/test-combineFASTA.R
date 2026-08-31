test_that("combineFASTA combines sequences from multiple FASTA files and saves to disk", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()

  f1 <- file.path(indir, "file1.fasta")
  f2 <- file.path(indir, "file2.fasta")
  fastadframe(fixture_seq_dframe()[1:2, ], file = f1)
  fastadframe(fixture_seq_dframe()[3, , drop = FALSE], file = f2)

  result <- combineFASTA(input_files = c(f1, f2), verbose = FALSE, dir = outdir)

  expect_s3_class(result, "combinedFASTA")
  # Note: total_sequences here reflects length(all_sequences), which is one
  # list element per INPUT FILE (each a DNAbin list), not the sum of
  # individual sequence records across files.
  expect_equal(result$summary$total_sequences, 2)
  expect_equal(result$summary$total_files, 2)
  expect_true(file.exists(result$output_path))
})

test_that("combineFASTA save = FALSE returns results without writing to disk", {
  indir <- withr::local_tempdir()

  f1 <- file.path(indir, "file1.fasta")
  fastadframe(fixture_seq_dframe(), file = f1)

  result <- combineFASTA(input_files = f1, save = FALSE, verbose = FALSE)

  expect_null(result$output_path)
  expect_equal(result$summary$output_file, "Not saved to disk")
})

test_that("combineFASTA errors when input_files is missing", {
  expect_error(combineFASTA(input_files = NULL), "required")
})

test_that("combineFASTA errors when input_files is not a character vector", {
  expect_error(combineFASTA(input_files = 123), "character vector")
})

test_that("combineFASTA errors when a specified file does not exist", {
  expect_error(combineFASTA(input_files = "does_not_exist.fasta"), "do not exist")
})

test_that("combineFASTA prints progress messages when verbose = TRUE", {
  indir <- withr::local_tempdir()
  outdir <- withr::local_tempdir()
  f1 <- file.path(indir, "file1.fasta")
  fastadframe(fixture_seq_dframe(), file = f1)

  expect_message(combineFASTA(input_files = f1, verbose = TRUE, dir = outdir),
                 "COMBINATION COMPLETE")
})

test_that("combineFASTA skips a FASTA file with no sequences and warns on unreadable files", {
  indir <- withr::local_tempdir()
  empty_file <- file.path(indir, "empty.fasta")
  writeLines(character(0), empty_file)
  f1 <- file.path(indir, "file1.fasta")
  fastadframe(fixture_seq_dframe(), file = f1)

  expect_message(result <- combineFASTA(input_files = c(empty_file, f1),
                                        verbose = TRUE, save = FALSE),
                 "No sequences found")
  expect_equal(result$summary$total_files, 1)
})
