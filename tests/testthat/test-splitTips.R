test_that("splitTips parses genus, species, voucher and genbank accession", {
  tips <- c("Genus_species_Cardoso1234_MN123456",
            "Genus_alter_Cardoso5678_MN654321")
  df <- splitTips(tiplabels = tips)

  expect_s3_class(df, "data.frame")
  expect_equal(df$genus, c("Genus", "Genus"))
  expect_equal(df$species, c("species", "alter"))
  expect_true("genbank" %in% names(df))
  expect_equal(df$genbank, c("MN123456", "MN654321"))
})

test_that("splitTips flags doubtful identifications (cf./aff.)", {
  tips <- c("Genus_cf_species_Cardoso1234_MN123456",
            "Genus_species_Cardoso5678_MN654321")
  df <- splitTips(tiplabels = tips)

  expect_true("doubtID" %in% names(df))
  expect_equal(df$doubtID[1], "cf")
  expect_true(is.na(df$doubtID[2]))
})

test_that("splitTips flags infraspecific epithets", {
  tips <- c("Genus_species_variety_Cardoso1234_MN123456",
            "Genus_alter_Cardoso5678_MN654321")
  df <- splitTips(tiplabels = tips)

  expect_true("infrasp" %in% names(df))
  expect_equal(df$infrasp[1], "variety")
  expect_true(is.na(df$infrasp[2]))
})

test_that("splitTips accepts space-delimited tip labels", {
  tips <- c("Genus species Cardoso1234 MN123456")
  df <- splitTips(tiplabels = tips)

  expect_equal(df$genus, "Genus")
  expect_equal(df$species, "species")
})

test_that("splitTips handles tip labels without a GenBank accession", {
  tips <- c("Genus_species_Cardoso1234", "Genus_alter_Cardoso5678")
  df <- splitTips(tiplabels = tips)

  expect_false("genbank" %in% names(df))
  expect_equal(df$genus, c("Genus", "Genus"))
  expect_equal(df$species, c("species", "alter"))
})

test_that("splitTips always includes the original tiplabels column (with underscores normalized to spaces)", {
  tips <- c("Genus_species_Cardoso1234_MN123456")
  df <- splitTips(tiplabels = tips)

  expect_equal(df$tiplabels, gsub("_", " ", tips))
})
