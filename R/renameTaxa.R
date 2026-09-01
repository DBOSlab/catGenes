#' Rename taxa across multiple DNA alignments
#'
#' @author Domingos Cardoso
#'
#' @description Reads every DNA alignment (NEXUS, FASTA or PHYLIP; any mix is
#' fine) found in a directory and renames the sequences whose names are found
#' in \code{lookup}, leaving any unmatched sequence name unchanged. This is
#' useful, for example, for replacing opaque voucher codes with a more
#' informative label (e.g. combining the taxon name with the original code),
#' after fishing out a subset of taxa with \code{\link{fishSeqs}}.
#'
#' @usage
#' renameTaxa(filepath = NULL,
#'            lookup = NULL,
#'            format = NULL,
#'            overwrite = FALSE,
#'            verbose = TRUE,
#'            save = TRUE,
#'            dir = "RESULTS_renameTaxa")
#'
#' @param filepath Path to the directory where the DNA alignments are stored.
#'
#' @param lookup Either a named character vector, where each name is the
#' original sequence name and each value is the new name to be assigned to
#' it, or a two-column data.frame/table with the original names in the first
#' column and the new names in the second. Any sequence name not found in
#' \code{lookup} is kept unchanged.
#'
#' @param format Define either "NEXUS", "FASTA" or "PHYLIP" for writing the
#' resulting renamed alignments. The default (\code{NULL}) keeps each output
#' file in the same format as its original input file. Ignored when
#' \code{overwrite = TRUE}, in which case the original format is always kept.
#'
#' @param overwrite Logical, if \code{TRUE} the renamed alignment overwrites
#' the original input file in place, keeping its original format. The
#' default (\code{FALSE}) saves the renamed alignments into \code{dir}
#' instead, leaving the original files untouched.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the analysis will not be printed in the console in full.
#'
#' @param save Logical, if \code{FALSE} the renamed alignments are returned
#' as an R object but not saved to disk (\code{overwrite} is then ignored).
#'
#' @param dir The path to the directory where the renamed alignments should
#' be saved when \code{overwrite = FALSE}. The default is to create a
#' directory named **RESULTS_renameTaxa** and save the files within a
#' subfolder named after the current date.
#'
#' @seealso \code{\link{fishSeqs}}
#'
#' @return An invisible named list of the renamed alignments, one
#' list-formatted alignment per input file (in the same style produced by
#' \code{\link[ape]{read.nexus.data}}).
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' lookup <- c(FC520 = "Moldenhawera_blanchetiana_FC520",
#'             FC459 = "Tachigali_costaricensis_FC459")
#'
#' renameTaxa(filepath = "RESULTS_fishSeqs/01Sep2026",
#'            lookup = lookup,
#'            dir = "RESULTS_renameTaxa")
#'
#' # Rename in place, overwriting the original files
#' renameTaxa(filepath = "RESULTS_fishSeqs/01Sep2026",
#'            lookup = lookup,
#'            overwrite = TRUE)
#'}
#'
#' @importFrom ape read.nexus.data read.FASTA
#'
#' @export
#'
renameTaxa <- function(filepath = NULL,
                       lookup = NULL,
                       format = NULL,
                       overwrite = FALSE,
                       verbose = TRUE,
                       save = TRUE,
                       dir = "RESULTS_renameTaxa") {

  if (is.null(filepath)) {
    stop("Please provide 'filepath', the directory where the DNA alignments are stored.")
  }
  if (is.null(lookup) || length(lookup) == 0) {
    stop("Please provide 'lookup', mapping original names to new names.")
  }

  input_files <- list.files(filepath)
  if (length(input_files) == 0) {
    stop(paste0("There is no DNA alignment in the directory.\n",
                "Make sure you have provided a correct filepath.\n\n"),
         "Find help also with:\n",
         "Domingos Cardoso (JBRJ; cardosobot@gmail.com)")
  }

  if (save && !overwrite) {
    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername, recursive = TRUE)
    }
  }

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

    renamed <- .rename_seqs(seqs, lookup)

    if (verbose) {
      message("  Renamed ", attr(renamed, "n_renamed"), " of ", length(renamed),
             " sequence name(s); ", attr(renamed, "n_unmatched"),
             " left unchanged.")
    }

    outname <- gsub("[.][^.]*$", "", inputfile)

    if (save) {
      if (overwrite) {
        out_format <- .detect_alignment_format(filepathfile)
        outfile <- filepathfile
      } else {
        out_format <- if (!is.null(format)) format else .detect_alignment_format(filepathfile)
        ext <- switch(out_format, NEXUS = ".nex", PHYLIP = ".phy", ".fasta")
        outfile <- paste0(foldername, "/", outname, ext)
      }

      if (out_format == "NEXUS") {
        nexusdframe(renamed, file = outfile)
      } else if (out_format == "PHYLIP") {
        phylipdframe(renamed, file = outfile)
      } else {
        fastadframe(renamed, file = outfile)
      }
    }

    attributes(renamed)[c("n_renamed", "n_unmatched")] <- NULL
    results[[outname]] <- renamed
  }

  return(invisible(results))
}


#-------------------------------------------------------------------------------
# Auxiliary function to rename the elements of a list-formatted alignment
# (as produced by .read_alignment_list()/ape::read.nexus.data()) using a
# lookup table. Any name not found in lookup is left unchanged. Reports how
# many names were matched/unmatched as attributes on the returned list.

.rename_seqs <- function(seqs, lookup) {

  if (is.data.frame(lookup)) {
    lookup <- stats::setNames(as.character(lookup[[2]]), as.character(lookup[[1]]))
  }

  if (is.null(names(lookup)) || any(!nzchar(names(lookup)))) {
    stop("'lookup'/'rename' must be a named vector (or two-column table) ",
        "mapping original names to new names.")
  }

  old_names <- names(seqs)
  matched <- old_names %in% names(lookup)

  new_names <- old_names
  new_names[matched] <- lookup[old_names[matched]]

  names(seqs) <- new_names

  attr(seqs, "n_renamed") <- sum(matched)
  attr(seqs, "n_unmatched") <- sum(!matched)

  return(seqs)
}
