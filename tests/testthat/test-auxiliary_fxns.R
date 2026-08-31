test_that(".namedlist captures argument names as list names", {
  gene1 <- 1:3
  gene2 <- 4:6
  out <- catGenes:::.namedlist(gene1, gene2)

  expect_named(out, c("gene1", "gene2"))
  expect_equal(out$gene1, gene1)
  expect_equal(out$gene2, gene2)
})

test_that(".shortaxlabels strips the trailing identifier keeping genus_species", {
  df <- data.frame(species = c("Genus_species_Voucher1", "Genus_alter_Voucher2"),
                   sequence = c("ACGT", "TTGG"),
                   stringsAsFactors = FALSE)
  out <- catGenes:::.shortaxlabels(df)

  expect_equal(out$species, c("Genus_species", "Genus_alter"))
})

test_that(".shortaxlabelsmult strips only the trailing GenBank-like identifier", {
  df <- data.frame(species = c("Genus_species_Coll1_MN123456", "Genus_alter_Coll2_MN654321"),
                   sequence = c("ACGT", "TTGG"),
                   stringsAsFactors = FALSE)
  out <- catGenes:::.shortaxlabelsmult(df)

  expect_equal(out$species, c("Genus_species_Coll1", "Genus_alter_Coll2"))
})

test_that(".intersectAll finds the common elements across multiple vectors", {
  expect_equal(sort(catGenes:::.intersectAll(list(c(1, 2, 3), c(2, 3, 4)))), c(2, 3))
  expect_equal(sort(catGenes:::.intersectAll(list(c(1, 2, 3), c(2, 3, 4), c(3, 4, 5)))), 3)
})

test_that(".intersectAll errors with fewer than two vectors", {
  expect_error(catGenes:::.intersectAll(list(c(1, 2, 3))))
})

test_that("equalnumb detects whether all values in a vector are equal", {
  expect_true(catGenes:::equalnumb(c(5, 5, 5)))
  expect_true(catGenes:::equalnumb(numeric(0)))
  expect_false(catGenes:::equalnumb(c(5, 6, 5)))
  expect_true(catGenes:::equalnumb(c(5, NA, 5)))
})

test_that(".seq_revcompl reverse-complements a DNA sequence", {
  expect_equal(catGenes:::.seq_revcompl("ACGT", strand = "complement"), "ACGT")
  expect_equal(catGenes:::.seq_revcompl("AACG", strand = "complement"), "CGTT")
  # Any strand other than "complement" skips the reversal but the base-by-base
  # complement is still applied unconditionally.
  expect_equal(catGenes:::.seq_revcompl("AACG", strand = "plus"), "TTGC")
})

test_that(".seq_revcompl handles RNA sequences (uses U instead of T)", {
  expect_equal(catGenes:::.seq_revcompl("AACG", strand = "complement"), "CGTT")
  expect_equal(catGenes:::.seq_revcompl("AACGU", strand = "complement"), "ACGUU")
})

test_that(".replace_terminal_gaps converts leading/trailing dashes to '?'", {
  df <- data.frame(species = c("sp1", "sp2"),
                   sequence = c("--ACGT--", "ACGT"),
                   stringsAsFactors = FALSE)
  out <- catGenes:::.replace_terminal_gaps(df)

  expect_equal(out$sequence[1], "??ACGT??")
  expect_equal(out$sequence[2], "ACGT")
})

test_that(".replace_terminal_gaps leaves internal gaps untouched", {
  df <- data.frame(species = "sp1",
                   sequence = "AC--GT",
                   stringsAsFactors = FALSE)
  out <- catGenes:::.replace_terminal_gaps(df)

  expect_equal(out$sequence, "AC--GT")
})

test_that(".tax_voucher_adjust cleans Species/Voucher columns of an inputdf", {
  df <- data.frame(Species = c("Genus species", "Genus  alter"),
                   Voucher = c("Cardoso 1234", NA),
                   ITS = c(" MN111111 ", "MN222222"),
                   stringsAsFactors = FALSE)
  out <- catGenes:::.tax_voucher_adjust(inputdf = df, gb.colnames = "ITS")

  expect_equal(out$Species, c("Genus_species", "Genus_alter"))
  expect_equal(out$Voucher, c("Cardoso1234", "Unvouchered"))
  expect_equal(out$ITS, c("MN111111", "MN222222"))
})

test_that(".tax_voucher_adjust also works with multiple gene columns", {
  df <- data.frame(Species = c("Genus species", "Genus alter"),
                   ITS = c(" MN111111 ", "MN222222"),
                   matK = c("MN333333 ", " MN444444"),
                   stringsAsFactors = FALSE)
  out <- catGenes:::.tax_voucher_adjust(inputdf = df, gb.colnames = c("ITS", "matK"))

  expect_equal(out$ITS, c("MN111111", "MN222222"))
  expect_equal(out$matK, c("MN333333", "MN444444"))
})

test_that(".tax_voucher_adjust cleans standalone taxon/voucher/genbank vectors", {
  out <- catGenes:::.tax_voucher_adjust(inputdf = NULL,
                                        taxon = "Genus species",
                                        voucher = "Cardoso 1234",
                                        genbank = "MN 111111")

  expect_equal(out[[1]], "Genus_species")
  expect_equal(out[[2]], "Cardoso1234")
  expect_equal(out[[3]], "MN111111")
})

test_that(".tax_voucher_adjust marks an all-unvouchered vector as NA", {
  out <- catGenes:::.tax_voucher_adjust(inputdf = NULL,
                                        taxon = NULL,
                                        voucher = c(NA, NA),
                                        genbank = c("MN111111", "MN222222"))

  expect_true(is.na(out[[2]]))
})
