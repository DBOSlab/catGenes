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

test_that("splitTips captures the epithet (not the rank marker) after var./subsp.", {
  tips <- c("Pickeringia_montana_var_montana_AY386863",
            "Genus_species_subsp_epithet_Cardoso1234_MN123456",
            "Genus species var. epithet2 MN000001",
            "Genus_alter_Cardoso5678_MN654321")
  df <- splitTips(tiplabels = tips)

  expect_equal(df$infrasp[1], "montana")
  expect_equal(df$infrasp[2], "epithet")
  expect_equal(df$infrasp[3], "epithet2")
  expect_true(is.na(df$infrasp[4]))

  expect_equal(df$genbank, c("AY386863", "MN123456", "MN000001", "MN654321"))
  expect_equal(df$voucher[2], "Cardoso1234")
  expect_true(is.na(df$voucher[1]))
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

test_that("splitTips errors when neither tiplabels nor file is provided", {
  expect_error(splitTips(), "provide either")
})

test_that("splitTips errors when file does not exist", {
  expect_error(splitTips(file = "no_such_file.nex"), "does not exist")
})

test_that("splitTips reads tip labels from a NEXUS file", {
  tips <- c("Genus_alpha_Cardoso1234_MN123456", "Genus_beta_Cardoso5678_MN654321")
  seqs <- data.frame(species = tips, sequence = c("ACGTACGTAC", "ACGAACGTAC"))
  nex <- withr::local_tempfile(fileext = ".nex")
  nexusdframe(seqs, file = nex)

  df <- splitTips(file = nex)

  expect_equal(sort(df$genus), c("Genus", "Genus"))
  expect_equal(sort(df$genbank), c("MN123456", "MN654321"))
})

test_that("splitTips reads tip labels from a FASTA file", {
  tips <- c("Genus_alpha_Cardoso1234_MN123456", "Genus_beta_Cardoso5678_MN654321")
  seqs <- data.frame(species = tips, sequence = c("ACGTACGTAC", "ACGAACGTAC"))
  fas <- withr::local_tempfile(fileext = ".fasta")
  fastadframe(seqs, file = fas)

  df <- splitTips(file = fas)

  expect_equal(sort(df$genbank), c("MN123456", "MN654321"))
})

test_that("splitTips reads tip labels from a PHYLIP file", {
  tips <- c("Genus_alpha_Cardoso1234_MN123456", "Genus_beta_Cardoso5678_MN654321")
  seqs <- data.frame(species = tips, sequence = c("ACGTACGTAC", "ACGAACGTAC"))
  phy <- withr::local_tempfile(fileext = ".phy")
  phylipdframe(seqs, file = phy)

  df <- splitTips(file = phy)

  expect_equal(sort(df$genbank), c("MN123456", "MN654321"))
})

test_that("splitTips ignores file when tiplabels is also provided", {
  tips <- c("Genus_species_Cardoso1234_MN123456")
  df <- splitTips(tiplabels = tips, file = "no_such_file.nex")

  expect_equal(df$genus, "Genus")
})
