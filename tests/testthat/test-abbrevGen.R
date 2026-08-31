test_that("abbrevGen abbreviates to 1 letter when all genera differ", {
  tips <- c("Alpha_species", "Beta_species", "Gamma_species")
  df <- abbrevGen(tiplabels = tips)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_equal(df$abbreviation, c("A.", "B.", "G."))
  expect_equal(df$abbrev_tiplabels, c("A. species", "B. species", "G. species"))
})

test_that("abbrevGen extends to 2 letters when first letters collide", {
  tips <- c("Alpha_species", "Amanita_species", "Beta_species")
  df <- abbrevGen(tiplabels = tips)

  # Because at least one pair of genera collides on the first letter, ALL
  # genera in the call are abbreviated to 2 letters, not just the colliding ones.
  expect_equal(unique(df$abbreviation[df$original_genus == "Alpha"]), "Al.")
  expect_equal(unique(df$abbreviation[df$original_genus == "Amanita"]), "Am.")
  expect_equal(unique(df$abbreviation[df$original_genus == "Beta"]), "Be.")
})

test_that("abbrevGen extends to 3 letters when first two letters collide", {
  tips <- c("Alaria_species", "Alamos_species", "Beta_species")
  df <- abbrevGen(tiplabels = tips)

  expect_equal(unique(df$abbreviation[df$original_genus == "Alaria"]), "Ala.")
  expect_equal(unique(df$abbreviation[df$original_genus == "Alamos"]), "Ala.")
})

test_that("abbrevGen accepts underscore- or space-delimited tip labels", {
  tips_space <- c("Alpha species voucher1", "Beta species voucher2")
  df <- abbrevGen(tiplabels = tips_space)
  expect_equal(df$original_genus, c("Alpha", "Beta"))
})

test_that("abbrevGen abbrevmult mode disambiguates progressively", {
  tips <- c("Alaria_species", "Alamos_species", "Beta_species")
  df <- abbrevGen(tiplabels = tips, abbrevfull = FALSE, abbrevmult = TRUE)

  expect_true(all(!is.na(df$abbreviation)))
  # Alaria and Alamos share the same first two letters, so must be disambiguated
  # from each other, while Beta (unique first letter) stays a single-letter abbreviation.
  expect_false(identical(unique(df$abbreviation[df$original_genus == "Alaria"]),
                         unique(df$abbreviation[df$original_genus == "Alamos"])))
  expect_equal(unique(df$abbreviation[df$original_genus == "Beta"]), "B.")
})
