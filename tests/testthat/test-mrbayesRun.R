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

test_that("mrbayesRun errors when nexus_file does not exist", {
  expect_error(mrbayesRun("does_not_exist.nex", mrbayes_dir = "."), "does not exist")
})

test_that("mrbayesRun errors when mrbayes_dir does not exist", {
  nex <- .make_fake_nexus(environment())
  expect_error(mrbayesRun(nex, mrbayes_dir = "no_such_dir_xyz"), "does not exist")
})

test_that("mrbayesRun errors when the mb executable is not found", {
  nex <- .make_fake_nexus(environment())
  empty_dir <- withr::local_tempdir()
  expect_error(mrbayesRun(nex, mrbayes_dir = empty_dir), "executable not found")
})

test_that("mrbayesRun errors on invalid copy_mode", {
  expect_error(mrbayesRun("x.nex", mrbayes_dir = ".", copy_mode = "bogus"),
              "copy_mode must be")
})

test_that("mrbayesRun runs via the base-R system2 fallback and copies the nexus file", {
  nex <- .make_fake_nexus(environment())
  mrbayes_dir <- .make_fake_mb_dir(environment())
  run_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(system2 = function(...) 0L, .package = "base")

  res <- mrbayesRun(nex,
                    mrbayes_dir = mrbayes_dir,
                    run_dir = run_dir,
                    use_processx = FALSE,
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
                    mrbayes_dir = mrbayes_dir,
                    run_dir = run_dir,
                    copy_mode = "move",
                    use_processx = FALSE,
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
    read_output_lines = function() character(0),
    read_error_lines = function() character(0),
    get_exit_status = function() 0L,
    get_pid = function() 1234L,
    signal = function(sig) invisible(TRUE),
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
                    mrbayes_dir = mrbayes_dir,
                    run_dir = run_dir,
                    use_processx = TRUE,
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
                    mrbayes_dir = mrbayes_dir,
                    run_dir = run_dir,
                    use_processx = TRUE,
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
