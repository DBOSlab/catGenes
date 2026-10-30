#' Run MrBayes from RStudio using an existing NEXUS file
#'
#' @author Domingos Cardoso
#'
#' @description
#' Run a Bayesian phylogenetic analysis in
#' [MrBayes](https://nbisweden.github.io/MrBayes/) directly from R/RStudio using a
#' pre-configured NEXUS file that already contains a complete \code{begin mrbayes;}
#' command block. The function copies (or moves) the NEXUS file to a dedicated run
#' directory, launches MrBayes, and records the program output in log files.
#'
#' When \pkg{processx} is installed and \code{use.processx = TRUE}, output can be
#' streamed live to the R console and the returned object includes a \code{stop()}
#' handle that can terminate the run programmatically. Programmatic stopping is
#' only useful when \code{background = TRUE}, because in foreground mode the call
#' blocks until MrBayes finishes.
#'
#' Without \pkg{processx}, the analysis can still be run and live output can be
#' printed, but stopping is interactive (e.g., Esc in RStudio).
#'
#' @details
#' Before starting, the NEXUS file is checked for non-ASCII characters (e.g.,
#' accented letters in collector names inside \code{[...]} comments). MrBayes cannot read them and stops with errors such as
#' \emph{"Reached end of file while in comment"}, so the function stops early and
#' reports the offending lines.
#'
#' Output capture:
#' - \code{stdout}: file \code{mrbayes_stdout.txt} in \code{run.dir}, where the
#'   MrBayes screen output is recorded.
#' - \code{stderr}: file \code{mrbayes_stderr.txt} in \code{run.dir}, where errors
#'   and warnings are recorded.
#'
#' With \pkg{processx}, MrBayes writes directly to these files, so they grow while
#' the analysis runs, also in background mode, whether or not \code{poll()} is
#' called. The function runs MrBayes in \code{run.dir}, so all MrBayes output files
#' (e.g., \code{.p}, \code{.t}, \code{.con.tre}, named after the \code{filename}
#' option of the \code{mcmc} command or, if absent, after the NEXUS file) are
#' written there.
#'
#' Stopping the run:
#' - With \pkg{processx} + \code{background = TRUE}: call \code{mrbayesStop(res)}
#'   or \code{res$stop()}. This interrupts MrBayes and confirms the stop, so the
#'   samples written so far are kept. Use \code{mrbayesStop(res, force = TRUE)} to
#'   kill the process immediately.
#' - With \pkg{processx} + \code{background = FALSE}: stop interactively (Esc) because the call blocks.
#' - Without \pkg{processx}: stop interactively (Esc in RStudio).
#'
#' Avoid placing \code{run.dir} inside folders synchronized by cloud services
#' (e.g., iCloud Drive, Dropbox), which may lock or offload files while MrBayes is
#' writing them during long runs.
#'
#' @usage
#' mrbayesRun(nexus.file,
#'            mrbayes.dir,
#'            run.dir = NULL,
#'            copy.mode = "copy",
#'            quiet = FALSE,
#'            live = TRUE,
#'            use.processx = TRUE,
#'            background = FALSE)
#'
#' @param nexus.file
#' Path to a NEXUS file (\code{.nex} or \code{.nexus}) that already contains a valid
#' MrBayes block (i.e., \code{begin mrbayes;} ... \code{end;}).
#'
#' @param mrbayes.dir
#' Path to the folder containing the MrBayes executable. On macOS/Linux this is
#' typically a folder containing \code{mb}; on Windows, \code{mb.exe}. If the
#' executable is not found directly in \code{mrbayes.dir}, the function also checks
#' \code{file.path(mrbayes.dir, "bin")}. Typical values for \code{mrbayes.dir}:
#' - macOS (Intel, Homebrew): \code{"/usr/local/bin"}
#' - macOS (Apple Silicon, Homebrew): \code{"/opt/homebrew/bin"}
#' - macOS/Linux (compiled from source): the \code{src} folder of the MrBayes
#'   source tree, e.g. \code{"~/MrBayes/src"}
#' - macOS/Linux (local install): \code{"~/.local/bin"} (use \code{path.expand("~/.local/bin")})
#' - Windows (typical): \code{"C:/Program Files/MrBayes"} or a folder containing \code{mb.exe}
#'
#' @param run.dir
#' Directory where the analysis will run. If \code{NULL}, a new folder named
#' \code{mrbayes_run_YYYYmmdd_HHMMSS} is created in \code{getwd()}.
#'
#' @param copy.mode
#' Either \code{"copy"} (default) or \code{"move"}. Determines whether \code{nexus.file}
#' is copied to \code{run.dir} or moved there before running MrBayes.
#'
#' @param quiet
#' Logical. If \code{TRUE}, suppresses most informational messages (default: \code{FALSE}).
#'
#' @param live
#' Logical. If \code{TRUE}, prints MrBayes output to the R console (default: \code{TRUE}).
#' In background mode, output is printed only when \code{poll()} is called.
#'
#' @param use.processx
#' Logical. If \code{TRUE} (default), uses \pkg{processx} when installed to provide
#' robust output streaming and a programmatic stop handle.
#'
#' @param background
#' Logical. If \code{TRUE}, starts MrBayes and returns immediately (non-blocking).
#' The returned object includes \code{poll()} to print new output on demand and
#' \code{stop()} to terminate the run. Requires \pkg{processx}. Default: \code{FALSE}.
#'
#' @return
#' A list with (at minimum) \code{exe_path}, \code{run_dir}, \code{nexus_path},
#' \code{stdout}, \code{stderr}, \code{exit_status}, and \code{error}. When
#' \pkg{processx} is used, the list additionally includes \code{proc} (a
#' \pkg{processx} process object) and \code{stop}. When \code{background = TRUE},
#' the list also includes \code{poll}, and \code{exit_status} is \code{NA} (use
#' \code{res$proc$is_alive()} and \code{res$proc$get_exit_status()} to follow the
#' run).
#'
#' @examples
#' \dontrun{
#' # Background mode (recommended for stop-anytime control)
#' res <- mrbayesRun("analysis.nex",
#'                   mrbayes.dir = "/usr/local/bin",
#'                   background = TRUE)
#' res$poll()             # print new output
#' res$proc$is_alive()    # still running?
#' mrbayesStop(res)       # graceful stop
#' mrbayesStop(res, TRUE) # force kill
#'
#' # Foreground mode (blocking): stop with Esc in RStudio
#' mrbayesRun("analysis.nex",
#'            mrbayes.dir = "/usr/local/bin")
#' }
#'
#' @seealso \code{\link{evomodelTest}}
#'
#' @importFrom processx process
#'
#' @export
#'

mrbayesRun <- function(nexus.file,
                       mrbayes.dir,
                       run.dir = NULL,
                       copy.mode = "copy",
                       quiet = FALSE,
                       live = TRUE,
                       use.processx = TRUE,
                       background = FALSE) {

  if (!copy.mode %in% c("copy", "move")) stop("copy.mode must be 'copy' or 'move'.")
  if (!is.logical(background) || length(background) != 1) stop("background must be TRUE/FALSE.")
  if (!is.logical(live) || length(live) != 1) stop("live must be TRUE/FALSE.")
  if (!is.logical(use.processx) || length(use.processx) != 1) stop("use.processx must be TRUE/FALSE.")

  if (is.null(nexus.file) || !nzchar(nexus.file)) stop("nexus.file is required.")
  if (!file.exists(nexus.file)) stop("nexus.file does not exist: ", nexus.file)
  nexus.file <- normalizePath(nexus.file, winslash = "/", mustWork = TRUE)

  # MrBayes cannot parse non-ASCII bytes (e.g. accented letters in comments)
  nex_bytes <- readBin(nexus.file, "raw", n = file.size(nexus.file))
  if (any(nex_bytes > as.raw(0x7F))) {
    nex_lines <- readLines(nexus.file, warn = FALSE, encoding = "bytes")
    bad <- which(vapply(nex_lines, function(l) any(charToRaw(l) > as.raw(0x7F)),
                        logical(1), USE.NAMES = FALSE))
    stop("nexus.file contains non-ASCII characters (e.g., accented letters), ",
         "which MrBayes cannot read.\n",
         "Lines: ", paste(utils::head(bad, 20), collapse = ", "),
         if (length(bad) > 20) ", ..." else "", "\n",
         "Replace them with plain ASCII (e.g., \u00e1 -> a) and run again.")
  }

  if (is.null(mrbayes.dir) || !nzchar(mrbayes.dir)) stop("mrbayes.dir is required.")
  if (!dir.exists(mrbayes.dir)) stop("mrbayes.dir does not exist: ", mrbayes.dir)
  mrbayes.dir <- normalizePath(mrbayes.dir, winslash = "/", mustWork = TRUE)

  exe_name <- if (identical(.Platform$OS.type, "windows")) "mb.exe" else "mb"

  exe_path <- file.path(mrbayes.dir, exe_name)
  if (!file.exists(exe_path)) {
    exe_alt <- file.path(mrbayes.dir, "bin", exe_name)
    if (!file.exists(exe_alt)) {
      stop("MrBayes executable not found.\nTried: ", exe_path, "\nAlso tried: ", exe_alt)
    }
    exe_path <- exe_alt
  }
  exe_path <- normalizePath(exe_path, winslash = "/", mustWork = TRUE)

  if (is.null(run.dir)) {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run.dir <- file.path(getwd(), paste0("mrbayes_run_", stamp))
  }
  dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)
  run.dir <- normalizePath(run.dir, winslash = "/", mustWork = TRUE)

  dest_nexus <- file.path(run.dir, basename(nexus.file))
  if (!identical(dest_nexus, nexus.file)) {
    if (copy.mode == "copy") {
      ok <- file.copy(nexus.file, dest_nexus, overwrite = TRUE)
      if (!ok) stop("Failed to copy nexus.file into run.dir.")
    } else {
      ok <- file.rename(nexus.file, dest_nexus)
      if (!ok) stop("Failed to move nexus.file into run.dir.")
    }
  }

  stdout_file <- file.path(run.dir, "mrbayes_stdout.txt")
  stderr_file <- file.path(run.dir, "mrbayes_stderr.txt")

  if (!quiet) {
    message("MrBayes executable: ", exe_path)
    message("Run directory: ", run.dir)
    message("Running NEXUS: ", dest_nexus)
  }

  want_processx <- isTRUE(use.processx) && requireNamespace("processx", quietly = TRUE)
  if (isTRUE(background) && !want_processx) stop("background=TRUE requires the processx package.")

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(run.dir)

  # ---- processx path (supports background + programmatic stop) ----
  if (want_processx) {
    # MrBayes writes straight to the log files, so its output never fills a pipe
    # (which would block the run) and the logs grow even if poll() is never called.
    # stdin stays open so that stop() can answer MrBayes's
    # "Do you really want to stop the run (y/n)?" prompt.
    proc <- processx::process$new(
      command = exe_path,
      args = c(basename(dest_nexus)),
      wd = run.dir,
      stdin = "|",
      stdout = stdout_file,
      stderr = stderr_file,
      cleanup = TRUE
    )

    stop_fun <- function(grace = TRUE) {
      if (!proc$is_alive()) return(invisible(FALSE))
      if (grace) {
        proc$interrupt()
        try(proc$write_input("y\n"), silent = TRUE)
        proc$wait(10000)
        if (proc$is_alive()) proc$kill()
      } else {
        proc$kill()
      }
      invisible(TRUE)
    }

    read_pos <- c(out = 0, err = 0)
    read_new <- function(path, key) {
      size <- file.size(path)
      if (is.na(size) || size <= read_pos[[key]]) return(character(0))
      con <- file(path, open = "rb")
      on.exit(close(con))
      seek(con, read_pos[[key]])
      txt <- rawToChar(readBin(con, "raw", n = size - read_pos[[key]]))
      read_pos[[key]] <<- size
      txt
    }

    poll_fun <- function() {
      out <- read_new(stdout_file, "out")
      if (length(out) && isTRUE(live)) cat(out)
      errl <- read_new(stderr_file, "err")
      if (length(errl) && isTRUE(live)) cat(errl, file = stderr())
      invisible(NULL)
    }

    if (!quiet) {
      msg <- if (background) {
        "Process started (background). Use res$poll() to follow; res$stop() to stop."
      } else {
        "Process started (foreground). Stop with Esc, or use background=TRUE for stop-anytime control."
      }
      message(msg, " (PID ", proc$get_pid(), ")")
    }

    if (isTRUE(background)) {
      return(list(
        exe_path = exe_path,
        run_dir = run.dir,
        nexus_path = dest_nexus,
        stdout = stdout_file,
        stderr = stderr_file,
        exit_status = NA_integer_,
        error = NULL,
        proc = proc,
        stop = stop_fun,
        poll = poll_fun
      ))
    }

    # Foreground (blocking) streaming
    exit_status <- NA_integer_
    err <- NULL
    tryCatch({
      while (proc$is_alive()) {
        poll_fun()
        Sys.sleep(0.2)
      }
      poll_fun()
      exit_status <- proc$get_exit_status()
    }, interrupt = function(e) {
      err <<- e
      try(proc$kill(), silent = TRUE)
      exit_status <<- NA_integer_
      if (!quiet) message("\nInterrupted: process killed.")
    }, error = function(e) {
      err <<- e
      try(proc$kill(), silent = TRUE)
      exit_status <<- NA_integer_
    })

    return(list(
      exe_path = exe_path,
      run_dir = run.dir,
      nexus_path = dest_nexus,
      stdout = stdout_file,
      stderr = stderr_file,
      exit_status = exit_status,
      error = err,
      proc = proc,
      stop = stop_fun
    ))
  }

  # ---- base R fallback (no process handle) ----
  if (!quiet && isTRUE(live)) {
    message("Tip: install.packages('processx') for background mode and programmatic stop().")
  }

  res <- tryCatch(
    system2(
      exe_path,
      args = c(basename(dest_nexus)),
      stdout = if (live) "" else stdout_file,
      stderr = if (live) "" else stderr_file
    ),
    error = function(e) structure(NA_integer_, error = e)
  )

  list(
    exe_path = exe_path,
    run_dir = run.dir,
    nexus_path = dest_nexus,
    stdout = stdout_file,
    stderr = stderr_file,
    exit_status = if (is.numeric(res)) res else NA_integer_,
    error = attr(res, "error")
  )
}

#' Stop a MrBayes run started by \code{mrbayesRun}
#'
#' @description
#' Convenience wrapper to stop a MrBayes run returned by \code{\link{mrbayesRun}}.
#' For programmatic stopping, the run must have been started with \pkg{processx}
#' (i.e., \code{use.processx = TRUE}) and, for practical "stop anytime" control,
#' \code{background = TRUE}.
#'
#' A graceful stop (\code{force = FALSE}) interrupts MrBayes and answers "yes" to
#' its \emph{"Do you really want to stop the run (y/n)?"} prompt, so the run ends
#' cleanly and the samples written so far are kept. If MrBayes has not exited
#' after 10 seconds, it is killed.
#'
#' @param res A result object returned by \code{\link{mrbayesRun}}.
#' @param force Logical. If \code{FALSE} (default), requests a graceful stop.
#' If \code{TRUE}, force-kills the process.
#'
#' @return Invisibly returns \code{TRUE} if a stop signal was sent, otherwise \code{FALSE}.
#'
#' @export
mrbayesStop <- function(res, force = FALSE) {
  if (is.null(res) || !is.list(res)) stop("res must be an object returned by mrbayesRun().")
  if (!is.logical(force) || length(force) != 1) stop("force must be TRUE/FALSE.")

  if (!"stop" %in% names(res) || !is.function(res$stop)) {
    stop(
      "No programmatic stop handle.\n",
      "Use processx and run: mrbayesRun(..., background = TRUE)\n",
      "Otherwise, stop interactively (Esc in RStudio)."
    )
  }

  invisible(res$stop(!isTRUE(force)))
}
