test_that(".arg_check_dir accepts character and strips trailing slash", {
  expect_equal(catGenes:::.arg_check_dir("some_dir/"), "some_dir")
  expect_equal(catGenes:::.arg_check_dir("some_dir"), "some_dir")
})

test_that(".arg_check_dir errors on non-character input", {
  expect_error(catGenes:::.arg_check_dir(123), "should be a character")
})

test_that(".arg_check_inputdf passes silently when no accessions are duplicated", {
  df <- data.frame(ITS = c("MN111111", "MN222222"),
                   matK = c("MN333333", "MN444444"),
                   stringsAsFactors = FALSE)
  expect_silent(catGenes:::.arg_check_inputdf(df, c("ITS", "matK")))
})

test_that(".arg_check_inputdf errors listing duplicated accessions", {
  df <- data.frame(ITS = c("MN111111", "MN111111"),
                   matK = c("MN333333", "MN444444"),
                   stringsAsFactors = FALSE)
  expect_error(catGenes:::.arg_check_inputdf(df, c("ITS", "matK")),
              "duplicated")
})

test_that(".arg_check_inputdf ignores NA when checking duplicates", {
  df <- data.frame(ITS = c(NA, NA, "MN111111"),
                   stringsAsFactors = FALSE)
  expect_silent(catGenes:::.arg_check_inputdf(df, "ITS"))
})
