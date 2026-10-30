#' Convert any format of DNA alignment into another format
#'
#' @author Domingos Cardoso
#'
#' @description Convert one or multiple DNA alignments in FASTA, NEXUS or PHYLIP
#' format into each other.
#'
#' @usage
#' convertAlign(filepath = NULL,
#'              format = NULL,
#'              rmfiles = FALSE,
#'              verbose = TRUE,
#'              dir = "RESULTS_convertAlign")
#'
#' @param filepath Path to the directory where the DNA alignments are stored.
#'
#' @param format Define either "NEXUS", "FASTA" or "PHYLIP" format for writing
#' the resulting newly converted files.
#'
#' @param rmfiles Logical, if \code{TRUE}, the original input file(s) are removed
#' from the directory and only the newly converted files are kept.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the analysis will not be printed in the console in full.
#'
#' @param dir The path to the directory where the newly converted files should
#' be saved. If no directory is given, then the files will saved in a subfolder
#' named after the current date under the folder named **RESULTS_convertAlign**.
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(GenBank_accessions)
#'
#' todaydate <- format(Sys.time(), "%d%b%Y")
#' folder_mined_seqs <- paste0("RESULTS_mineSeq/", todaydate)
#' folder_aligned_seqs <- paste0("RESULTS_alignSeqs/", todaydate)
#'
#' mineSeq(inputdf = GenBank_accessions,
#'         gb.colnames = c("ITS", "matK"),
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         dir = "RESULTS_mineSeq",
#'         filename = "GenBanK_seqs")
#'
#' alignSeqs(filepath = folder_mined_seqs,
#'           method = "ClustalW",
#'           gapOpening = "default",
#'           format = "NEXUS",
#'           verbose = TRUE,
#'           dir = "RESULTS_alignSeqs")
#'
#' convertAlign(filepath = folder_aligned_seqs,
#'           format = "PHYLIP",
#'           rmfiles = FALSE,
#'           verbose = TRUE,
#'           dir = "RESULTS_convertAlign")
#'}
#'
#' @importFrom ape read.nexus.data
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
convertAlign <- function(filepath = NULL,
                         format = NULL,
                         rmfiles = FALSE,
                         verbose = TRUE,
                         dir = "RESULTS_convertAlign") {

  if (is.null(filepath) || !dir.exists(filepath)) {
    stop("'filepath' must be the path to an existing directory of alignments.")
  }
  if (is.null(format) || !toupper(format) %in% c("NEXUS", "FASTA", "PHYLIP")) {
    stop("'format' must be one of 'NEXUS', 'FASTA' or 'PHYLIP'.")
  }
  format <- toupper(format)

  # Only regular files (skip sub-directories such as a previous output folder)
  inputfiles <- list.files(filepath)
  inputfiles <- inputfiles[!dir.exists(file.path(filepath, inputfiles))]
  # Drop only the extension so names like "gene.v2.fasta" stay distinct
  inputfiles_names <- tools::file_path_sans_ext(inputfiles)

  # Create folder to save the converted DNA sequences
  foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
  if (!dir.exists(foldername)) {
    dir.create(foldername, recursive = TRUE)
  }

  n_converted <- 0
  for (i in seq_along(inputfiles)) {

    infile <- file.path(filepath, inputfiles[i])
    # Non-text files (e.g. .xlsx, .RData) in the folder are skipped, not parsed
    temp <- tryCatch(suppressWarnings(readLines(infile, warn = FALSE)),
                     error = function(e) NULL)
    if (is.null(temp) || !all(validUTF8(temp))) {
      message("The input file ", inputfiles[i],
              " is not a text alignment and was skipped.\n")
      next
    }
    temp <- sub("^\ufeff", "", temp)
    temp <- temp[nzchar(trimws(temp))]
    if (length(temp) == 0) {
      message("The input file ", inputfiles[i], " is empty and was skipped.\n")
      next
    }
    firstchar <- substr(trimws(temp[1]), 1, 1)
    if (!firstchar %in% c("#", ">") &&
        !grepl("^\\s*[0-9]+\\s+[0-9]+", temp[1])) {
      message("The input file ", inputfiles[i],
              " is not a NEXUS, FASTA or PHYLIP alignment and was skipped.\n")
      next
    }

    if (firstchar == "#") {
      if (format == "NEXUS") {
        message(paste0("The input file ", inputfiles[i],
                       "\nis already in NEXUS format!\n"),
                "This file was not converted.\n\n")
        next
      }
      # Convert a NEXUS file ####
      filetoconvert <- ape::read.nexus.data(infile)
      for (j in seq_along(filetoconvert)) {
        filetoconvert[[j]] <- paste(filetoconvert[[j]], collapse = "")
        filetoconvert[[j]] <- toupper(filetoconvert[[j]])
      }
      spp <- stats::setNames(data.frame(names(filetoconvert)), "species")
      seqs <- stats::setNames(data.frame(unlist(filetoconvert, use.names = FALSE)), "sequence")
      convertedfile <- cbind(spp, seqs)

    } else if (firstchar == ">") {
      if (format == "FASTA") {
        message(paste0("The input file ", inputfiles[i],
                       "\nis already in FASTA format!\n"),
                "This file was not converted.\n\n")
        next
      }
      # Convert a FASTA file (sequences may be wrapped over several lines) ####
      tf <- grepl("^>", trimws(temp))
      recid <- cumsum(tf)
      convertedfile <- data.frame(
        species = trimws(sub("^>", "", trimws(temp[tf]))),
        sequence = vapply(split(trimws(temp[!tf]), recid[!tf]),
                          paste, character(1), collapse = ""),
        row.names = NULL)

    } else {
      if (format == "PHYLIP") {
        message(paste0("The input file ", inputfiles[i],
                       "\nis already in PHYLIP format!\n"),
                "This file was not converted.\n\n")
        next
      }
      # Convert a PHYLIP file ####
      filetoconvert <- temp[-1]
      convertedfile <- data.frame(species = gsub("\\s.*", "", filetoconvert),
                                  sequence = gsub(".*\\s", "", filetoconvert))
    }

    ext <- c(NEXUS = ".nex", PHYLIP = ".phy", FASTA = ".fasta")[[format]]
    outfile <- paste0(foldername, "/", inputfiles_names[i], ext)

    if (format == "NEXUS") {
      nexusdframe(convertedfile, outfile)
    } else if (format == "PHYLIP") {
      phylipdframe(convertedfile, outfile)
    } else if (format == "FASTA") {
      fastadframe(convertedfile, outfile)
    }

    n_converted <- n_converted + 1
    # Remove original file from the directory
    if (rmfiles) {
      unlink(infile)
    }
  }

  if (n_converted == 0) {
    warning("No alignment was converted, so nothing was saved in ", foldername,
            ". Input files already in ", format, " format are skipped.")
  } else if (verbose) {
    message(n_converted, " alignment(s) converted to ", format,
            " and saved in: ", foldername)
  }

  invisible(foldername)
}
