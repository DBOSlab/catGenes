#' Fish specific taxa out of multiple sequence alignments
#'
#' @author Domingos Cardoso
#'
#' @description Reads every DNA alignment (NEXUS, FASTA or PHYLIP; any mix is
#' fine) found in a directory, keeps only the sequences whose names match
#' \code{taxa}, and discards every other sequence. Any alignment site
#' (column) left entirely composed of gaps ("-") or missing data ("?") in the
#' retained sequences is then removed, so each output file is still a valid,
#' trimmed alignment. This is useful, for example, for pulling a focal
#' clade's sequences out of a large multi-locus set of per-gene alignments to
#' run a smaller, focused analysis.
#'
#' @usage
#' fishSeqs(filepath = NULL,
#'          taxa = NULL,
#'          rename = NULL,
#'          format = NULL,
#'          verbose = TRUE,
#'          save = TRUE,
#'          dir = "RESULTS_fishSeqs")
#'
#' @param filepath Path to the directory where the DNA alignments are stored.
#'
#' @param taxa A vector of one or more patterns to be matched (case-insensitive)
#' against each alignment's sequence names, so as to keep any sequence whose
#' name contains at least one of the given patterns. Providing a genus name,
#' for instance, keeps every sequence of that genus whenever tip labels
#' follow the usual \code{Genus_species_Voucher} convention.
#'
#' @param rename An optional lookup used to relabel the fished sequences,
#' applied after taxa are kept and gap-only columns removed. Either a named
#' character vector, where each name is the original sequence name and each
#' value is the new name, or a two-column data.frame/table with the original
#' names in the first column and the new names in the second (see
#' \code{\link{renameTaxa}}). Any sequence name not found in \code{rename} is
#' kept unchanged. The default (\code{NULL}) does not rename anything.
#'
#' @param format Define either "NEXUS", "FASTA" or "PHYLIP" for writing the
#' resulting fished alignments. The default (\code{NULL}) keeps each output
#' file in the same format as its original input file.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the analysis will not be printed in the console in full.
#'
#' @param save Logical, if \code{FALSE} the fished alignments are returned
#' as an R object but not saved to disk.
#'
#' @param dir The path to the directory where the fished alignments should
#' be saved. The default is to create a directory named **RESULTS_fishSeqs**
#' and save the files within a subfolder named after the current date.
#'
#' @return An invisible named list of the fished alignments, one
#' list-formatted alignment per input file (in the same style produced by
#' \code{\link[ape]{read.nexus.data}}). Files with none of the requested taxa
#' are skipped.
#'
#' @seealso \code{\link{renameTaxa}}
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' fishSeqs(filepath = "loci_alignments",
#'          taxa = c("FC520", "FC459"),
#'          rename = c(FC520 = "Moldenhawera_blanchetiana_FC520",
#'                     FC459 = "Tachigali_costaricensis_FC459"),
#'          dir = "RESULTS_fishSeqs")
#'}
#'
#' @importFrom ape read.nexus.data read.FASTA
#'
#' @export
#'
fishSeqs <- function(filepath = NULL,
                     taxa = NULL,
                     rename = NULL,
                     format = NULL,
                     verbose = TRUE,
                     save = TRUE,
                     dir = "RESULTS_fishSeqs") {

  if (is.null(filepath)) {
    stop("Please provide 'filepath', the directory where the DNA alignments are stored.")
  }
  if (is.null(taxa) || length(taxa) == 0) {
    stop("Please provide 'taxa', a vector of taxon name patterns to keep.")
  }

  input_files <- list.files(filepath)
  if (length(input_files) == 0) {
    stop(paste0("There is no DNA alignment in the directory.\n",
                "Make sure you have provided a correct filepath.\n\n"),
         "Find help also with:\n",
         "Domingos Cardoso (JBRJ; cardosobot@gmail.com)")
  }

  if (save) {
    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername, recursive = TRUE)
    }
  }

  taxa_pattern <- paste(taxa, collapse = "|")

  results <- list()

  for (i in seq_along(input_files)) {
    inputfile <- input_files[i]
    filepathfile <- paste0(filepath, "/", inputfile)

    if (verbose) {
      message("\nProcessing: ", inputfile, " (", i, "/", length(input_files), ")")
    }

    seqs <- tryCatch(.read_alignment_list(filepathfile), error = function(e) {
      warning("Could not read '", inputfile, "': ", e$message, call. = FALSE)
      NULL
    })

    if (is.null(seqs) || length(seqs) == 0) {
      if (verbose) message("  Skipped: could not be read as an alignment.")
      next
    }

    keep <- grepl(taxa_pattern, names(seqs), ignore.case = TRUE)

    if (!any(keep)) {
      if (verbose) message("  No matching taxa found in this file. Skipped.")
      next
    }

    fished <- seqs[keep]

    n_sites_before <- unique(lengths(fished))
    if (length(n_sites_before) > 1) {
      warning("Sequences in '", inputfile, "' are not all the same length; ",
             "skipping this file (is it really an alignment?).", call. = FALSE)
      next
    }

    fished <- .trim_allgap_columns(fished)

    if (is.null(fished)) {
      if (verbose) {
        message("  Skipped: no site with actual data remained for the matched taxa.")
      }
      next
    }

    if (verbose) {
      message("  Kept ", length(fished), " of ", length(seqs), " sequence(s); ",
             n_sites_before, " -> ", length(fished[[1]]), " site(s) after ",
             "removing gap-only columns.")
    }

    if (!is.null(rename)) {
      fished <- .rename_seqs(fished, rename)
      if (verbose) {
        message("  Renamed ", attr(fished, "n_renamed"), " of ", length(fished),
               " sequence name(s).")
      }
      attributes(fished)[c("n_renamed", "n_unmatched")] <- NULL
    }

    outname <- gsub("[.][^.]*$", "", inputfile)
    out_format <- if (!is.null(format)) format else .detect_alignment_format(filepathfile)

    if (save) {
      if (out_format == "NEXUS") {
        nexusdframe(fished, file = paste0(foldername, "/", outname, ".nex"))
      } else if (out_format == "PHYLIP") {
        phylipdframe(fished, file = paste0(foldername, "/", outname, ".phy"))
      } else {
        fastadframe(fished, file = paste0(foldername, "/", outname, ".fasta"))
      }
    }

    results[[outname]] <- fished
  }

  if (length(results) == 0) {
    warning("No alignment yielded any of the requested taxa.", call. = FALSE)
  }

  return(invisible(results))
}


#-------------------------------------------------------------------------------
# Auxiliary function to read a NEXUS, FASTA or PHYLIP alignment into the
# list-formatted representation produced by ape::read.nexus.data(), i.e. a
# named list where each element is a character vector of single-site
# characters. Format is sniffed from the first non-empty line, following the
# same convention used by convertAlign().

.read_alignment_list <- function(file) {

  temp <- readLines(file, warn = FALSE)
  temp <- temp[nzchar(trimws(temp))]

  first_char <- substr(trimws(temp[1]), 1, 1)

  if (first_char == "#") {
    # NEXUS
    seqs <- ape::read.nexus.data(file)

  } else if (first_char == ">") {
    # FASTA
    seqs <- as.character(ape::read.FASTA(file))

  } else {
    # PHYLIP: first line holds the ntax/nchar dimensions
    temp <- temp[-1]
    sp <- gsub("\\s.*", "", trimws(temp))
    sq <- gsub("^\\S+\\s+", "", trimws(temp))
    seqs <- stats::setNames(lapply(sq, function(s) strsplit(s, "")[[1]]), sp)
  }

  # Normalize case across formats: ape::read.FASTA() returns lowercase bases,
  # while NEXUS/PHYLIP text is used as-is, which would otherwise make casing
  # depend on the input format.
  seqs <- lapply(seqs, toupper)

  return(seqs)
}


#-------------------------------------------------------------------------------
# Auxiliary function to sniff the format of an alignment file, following the
# same convention used by convertAlign().

.detect_alignment_format <- function(file) {

  temp <- readLines(file, warn = FALSE)
  temp <- temp[nzchar(trimws(temp))]
  first_char <- substr(trimws(temp[1]), 1, 1)

  if (first_char == "#") {
    return("NEXUS")
  } else if (first_char == ">") {
    return("FASTA")
  } else {
    return("PHYLIP")
  }
}


#-------------------------------------------------------------------------------
# Auxiliary function to remove alignment sites (columns) that are entirely
# gaps ("-") or missing data ("?") across the provided sequences. Returns
# NULL if no site with actual data remains.

.trim_allgap_columns <- function(seqs) {

  mat <- do.call(rbind, seqs)

  allgap <- apply(mat, 2, function(col) all(grepl("^[-?]$", col)))

  if (any(allgap)) {
    mat <- mat[, !allgap, drop = FALSE]
  }

  if (ncol(mat) == 0) {
    return(NULL)
  }

  out <- lapply(seq_len(nrow(mat)), function(i) unname(mat[i, ]))
  names(out) <- names(seqs)

  return(out)
}
