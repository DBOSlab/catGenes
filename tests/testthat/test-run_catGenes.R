test_that("run_catGenes locates the bundled Shiny app directory", {
  appDir <- system.file("catGeneshiny", package = "catGenes")
  expect_true(nzchar(appDir))
  expect_true(dir.exists(appDir))
})

test_that("run_catGenes launches shiny::runApp with the bundled app directory", {
  captured <- list()
  testthat::local_mocked_bindings(
    runApp = function(appDir, launch.browser = TRUE, port = NULL) {
      captured$appDir <<- appDir
      captured$launch.browser <<- launch.browser
      captured$port <<- port
      invisible("fake-app-result")
    },
    .package = "shiny")

  result <- run_catGenes(launch.browser = FALSE, port = 1234)

  expect_equal(basename(captured$appDir), "catGeneshiny")
  expect_false(captured$launch.browser)
  expect_equal(captured$port, 1234)
  expect_equal(result, "fake-app-result")
})
