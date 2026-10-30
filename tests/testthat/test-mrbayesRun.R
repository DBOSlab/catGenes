# These tests exercise mrbayesRun()'s base-R fallback path (use_processx = FALSE)
# with a fake "mb" executable and a mocked system2(), so no real MrBayes binary
# or process is ever invoked.

.make_fake_mb_dir <- function(envir = parent.frame()) {
  mrbayes_dir <- withr::local_tempdir(.local_envir = envir)
  exe_name <- if (identical(.Platform$OS.type, "windows")) "mb.exe" else "mb"
  file.create(file.path(mrbayes_dir, exe_name))
  mrbayes_dir
}

.make_fake_nexus <- function(envir = parent.frame()) {
  nex <- withr::local_tempfile(fileext = ".nex", .local_envir = envir)
  writeLines(c("#NEXUS", "begin mrbayes;", "end;"), nex)
  nex
}

test_that("mrbayesRun errors when nexus.file does not exist", {
  expect_error(mrbayesRun("does_not_exist.nex",
                          mrbayes.dir = "."), "does not exist")
})

test_that("mrbayesRun errors when mrbayes.dir does not exist", {
  nex <- .make_fake_nexus(environment())
  expect_error(mrbayesRun(nex,
                          mrbayes.dir = "no_such_dir_xyz"), "does not exist")
})

test_that("mrbayesRun errors when the mb executable is not found", {
  nex <- .make_fake_nexus(environment())
  empty_dir <- withr::local_tempdir()
  expect_error(mrbayesRun(nex,
                          mrbayes.dir = empty_dir), "executable not found")
})

test_that("mrbayesRun errors on invalid copy.mode", {
  expect_error(
    mrbayesRun(
      "x.nex",
      mrbayes.dir = ".",
      copy.mode = "bogus"
    ),
    "copy\\.mode must be 'copy' or 'move'\\."
  )
})

test_that("mrbayesRun runs via the base-R system2 fallback and copies the nexus file", {
  nex <- .make_fake_nexus(environment())
  mrbayes_dir <- .make_fake_mb_dir(environment())
  run_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(system2 = function(...) 0L, .package = "base")

  res <- mrbayesRun(nex,
                    mrbayes.dir = mrbayes_dir,
                    run.dir = run_dir,
                    use.processx = FALSE,
                    quiet = TRUE,
                    live = FALSE)

  expect_equal(res$exit_status, 0L)
  expect_true(file.exists(res$nexus_path))
  expect_true(file.exists(file.path(res$run_dir, basename(nex))))
})

test_that("mrbayesRun copy_mode = 'move' removes the source nexus file", {
  nex <- .make_fake_nexus(environment())
  mrbayes_dir <- .make_fake_mb_dir(environment())
  run_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(system2 = function(...) 0L, .package = "base")

  res <- mrbayesRun(nex,
                    mrbayes.dir = mrbayes_dir,
                    run.dir = run_dir,
                    copy.mode = "move",
                    use.processx = FALSE,
                    quiet = TRUE,
                    live = FALSE)

  expect_false(file.exists(nex))
  expect_true(file.exists(res$nexus_path))
})

test_that("mrbayesStop errors when res has no programmatic stop handle", {
  expect_error(mrbayesStop(list(exit_status = 0L)), "No programmatic stop handle")
})

test_that("mrbayesStop errors when res is not a list", {
  expect_error(mrbayesStop("not-a-list"), "must be an object returned")
})

.fake_process_instance <- function() {
  calls <- 0
  list(
    is_alive = function() { calls <<- calls + 1; calls < 2 },
    get_exit_status = function() 0L,
    get_pid = function() 1234L,
    signal = function(sig) invisible(TRUE),
    interrupt = function() invisible(TRUE),
    write_input = function(str) invisible(TRUE),
    wait = function(timeout) invisible(NULL),
    kill = function() invisible(TRUE)
  )
}
.fake_process_generator <- list(new = function(...) .fake_process_instance())

test_that("mrbayesRun runs via processx in the foreground and returns a stop handle", {
  nex <- .make_fake_nexus(environment())
  mrbayes_dir <- .make_fake_mb_dir(environment())
  run_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(process = .fake_process_generator, .package = "processx")

  res <- mrbayesRun(nex,
                    mrbayes.dir = mrbayes_dir,
                    run.dir = run_dir,
                    use.processx = TRUE,
                    quiet = TRUE,
                    live = FALSE,
                    background = FALSE)

  expect_equal(res$exit_status, 0L)
  expect_true(is.function(res$stop))
  expect_null(res$poll)
})

test_that("mrbayesRun background = TRUE returns immediately with poll/stop handles", {
  nex <- .make_fake_nexus(environment())
  mrbayes_dir <- .make_fake_mb_dir(environment())
  run_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(process = .fake_process_generator, .package = "processx")

  res <- mrbayesRun(nex,
                    mrbayes.dir = mrbayes_dir,
                    run.dir = run_dir,
                    use.processx = TRUE,
                    quiet = TRUE,
                    live = FALSE,
                    background = TRUE)

  expect_true(is.na(res$exit_status))
  expect_true(is.function(res$poll))
  expect_true(is.function(res$stop))

  res$poll()
  expect_true(isTRUE(res$stop()))
})

test_that("mrbayesStop calls the stop handle with the expected grace flag", {
  called_with <- NULL
  fake_res <- list(stop = function(grace) {
    called_with <<- grace
    invisible(TRUE)
  })

  mrbayesStop(fake_res, force = FALSE)
  expect_true(called_with)

  mrbayesStop(fake_res, force = TRUE)
  expect_false(called_with)
})

test_that("mrbayesRun errors on non-ASCII characters and reports the lines", {
  nex <- withr::local_tempfile(fileext = ".nex")
  writeLines(c("#NEXUS", "[Gonz\u00e1lez 123]", "begin mrbayes;", "end;"), nex, useBytes = TRUE)
  mrbayes_dir <- .make_fake_mb_dir(environment())

  expect_error(mrbayesRun(nex, mrbayes.dir = mrbayes_dir, quiet = TRUE),
               "non-ASCII.*\nLines: 2\n")
})

# The tests below start a real child process: a small shell script standing in
# for "mb", so they run without MrBayes but exercise the actual processx wiring.
.make_script_mb_dir <- function(script, envir = parent.frame()) {
  mrbayes_dir <- withr::local_tempdir(.local_envir = envir)
  exe <- file.path(mrbayes_dir, "mb")
  writeLines(c("#!/bin/sh", script), exe)
  Sys.chmod(exe, "755")
  mrbayes_dir
}

test_that("background runs are not blocked by large output and poll() works", {
  skip_on_os("windows")
  skip_if_not_installed("processx")
  nex <- .make_fake_nexus(environment())
  run_dir <- withr::local_tempdir()
  # ~400 KB of output: far more than an OS pipe buffer (~64 KB)
  mrbayes_dir <- .make_script_mb_dir(
    "i=0; while [ $i -lt 8000 ]; do echo \"generation $i -- [-34564.957] (-34636.456) ..\"; i=$((i+1)); done"
  )

  res <- mrbayesRun(nex, mrbayes.dir = mrbayes_dir, run.dir = run_dir,
                    quiet = TRUE, live = FALSE, background = TRUE)
  res$proc$wait(20000)

  expect_false(res$proc$is_alive())
  expect_equal(res$proc$get_exit_status(), 0L)
  expect_length(readLines(res$stdout), 8000)
  expect_no_error(res$poll())
})

test_that("foreground runs record output in the stdout log", {
  skip_on_os("windows")
  skip_if_not_installed("processx")
  nex <- .make_fake_nexus(environment())
  run_dir <- withr::local_tempdir()
  mrbayes_dir <- .make_script_mb_dir("echo hello from mb; echo oops >&2")

  res <- mrbayesRun(nex, mrbayes.dir = mrbayes_dir, run.dir = run_dir,
                    quiet = TRUE, live = FALSE)

  expect_equal(res$exit_status, 0L)
  expect_equal(readLines(res$stdout), "hello from mb")
  expect_equal(readLines(res$stderr), "oops")
})

test_that("graceful stop interrupts the run and answers MrBayes's prompt", {
  skip_on_os("windows")
  skip_if_not_installed("processx")
  nex <- .make_fake_nexus(environment())
  run_dir <- withr::local_tempdir()
  # Mimics MrBayes: on SIGINT ask for confirmation on stdin, stop only on "y"
  mrbayes_dir <- .make_script_mb_dir(c(
    "trap 'echo \"Do you really want to stop the run (y/n)?\"; read a; [ \"$a\" = y ] && { echo stopped; exit 0; }' INT",
    "while true; do sleep 0.1; done"
  ))

  res <- mrbayesRun(nex, mrbayes.dir = mrbayes_dir, run.dir = run_dir,
                    quiet = TRUE, live = FALSE, background = TRUE)
  Sys.sleep(0.5)

  expect_true(mrbayesStop(res))
  expect_false(res$proc$is_alive())
  expect_equal(res$proc$get_exit_status(), 0L)
  expect_true("stopped" %in% readLines(res$stdout))
  expect_false(mrbayesStop(res))
})
