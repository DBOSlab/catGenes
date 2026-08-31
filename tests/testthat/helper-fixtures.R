# Shared small synthetic fixtures used across multiple test files.
# Mimic the structure of ape::read.nexus.data() output: a named list where
# each element is a character vector of single-base characters, and names
# follow the "Genus_species_Voucher" convention used throughout catGenes.

.make_gene <- function(names_seqs) {
  # names_seqs: named character vector, name = taxon label, value = sequence string
  out <- lapply(names_seqs, function(s) toupper(strsplit(s, "")[[1]]))
  names(out) <- names(names_seqs)
  out
}

# Two small, non-duplicated-species gene alignments for catfullGenes()
fixture_geneA_full <- function() {
  .make_gene(c(
    Genus_alpha_V1 = "ACGTACGTAC",
    Genus_beta_V2  = "ACGAACGTAC",
    Genus_gamma_V3 = "ACGTTCGTAC",
    Genus_delta_V4 = "ACGTACCTAC"
  ))
}

fixture_geneB_full <- function() {
  .make_gene(c(
    Genus_alpha_V1 = "TTGGCCAATT",
    Genus_beta_V2  = "TTGGCCAAAT",
    Genus_gamma_V3 = "TTGGCCTATT",
    Genus_delta_V4 = "TTGCCCAATT"
  ))
}

# Gene C only shares alpha/beta/gamma with A/B (delta missing) - useful for
# missdata = FALSE / outgroup tests.
fixture_geneC_partial <- function() {
  .make_gene(c(
    Genus_alpha_V1 = "GGGGCCCCAA",
    Genus_beta_V2  = "GGGACCCCAA",
    Genus_gamma_V3 = "GGGGCACCAA"
  ))
}

# Multi-accession gene alignments (duplicated species, distinct vouchers) for
# catmultGenes()/dropSeq() tests. Two accessions of Genus_alpha with differing
# amounts of missing data ("?"), one accession each of beta and gamma.
fixture_geneA_mult <- function() {
  .make_gene(c(
    Genus_alpha_Coll1 = "ACGTACGTAC",
    Genus_alpha_Coll2 = "ACGTACGT??",
    Genus_beta_Coll3  = "ACGAACGTAC",
    Genus_gamma_Coll4 = "ACGTTCGTAC"
  ))
}

fixture_geneB_mult <- function() {
  .make_gene(c(
    Genus_alpha_Coll1 = "TTGGCCAATT",
    Genus_alpha_Coll2 = "TTGGCCAA??",
    Genus_beta_Coll3  = "TTGGCCAAAT",
    Genus_gamma_Coll4 = "TTGGCCTATT"
  ))
}

# A ready two-column species/sequence data.frame, as consumed directly by
# nexusdframe()/phylipdframe()/fastadframe()/convertAlign() writers.
fixture_seq_dframe <- function() {
  data.frame(
    species = c("Genus_alpha", "Genus_beta", "Genus_gamma"),
    sequence = c("ACGTACGTAC", "ACGAACGTAC", "ACGTTCGTAC"),
    stringsAsFactors = FALSE
  )
}
