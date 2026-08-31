test_that("evomodelCmds groups identical models across partitions", {
  models <- c("GTR+G(4)+I", "HKY+G(4)", "GTR+G(4)+I")
  out <- evomodelCmds(model_vector = models)

  expect_named(out, c("grouped_commands", "all_commands", "partition_map", "raw_data"))
  # Partitions 1 and 3 share the same model, so should be grouped together
  expect_equal(nrow(out$grouped_commands), 2)
  expect_true("1,3" %in% out$grouped_commands$partitions)
  expect_true("2" %in% out$grouped_commands$partitions)
})

test_that("evomodelCmds produces valid lset/prset command syntax", {
  out <- evomodelCmds(model_vector = c("GTR+G", "HKY"))

  expect_true(any(grepl("^lset applyto=\\(1\\) nst=6 rates=gamma", out$all_commands)))
  expect_true(any(grepl("^lset applyto=\\(2\\) nst=2 rates=equal", out$all_commands)))
  expect_true(any(grepl("^prset applyto=\\(1\\)", out$all_commands)))
})

test_that("evomodelCmds partition_map reflects the original model per position", {
  out <- evomodelCmds(model_vector = c("JC", "GTR+I"))

  expect_equal(unname(out$partition_map), c("JC", "GTR+I"))
  expect_equal(names(out$partition_map), c("1", "2"))
})

test_that("evomodelCmds defaults an unrecognized model to GTR with a warning", {
  expect_warning(out <- evomodelCmds(model_vector = "NOTAMODEL"), "not recognized")
  expect_equal(out$raw_data$nst, 6)
})

test_that("evomodelCmds errors on empty or non-character input", {
  expect_error(evomodelCmds(model_vector = character(0)), "non-empty character")
  expect_error(evomodelCmds(model_vector = 123), "non-empty character")
})

test_that("evomodelCmds handles invariant-sites-only models", {
  out <- evomodelCmds(model_vector = "HKY+I")
  expect_equal(out$raw_data$rates, "propinv")
})
