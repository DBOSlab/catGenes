test_that("evomodelCmds groups identical models across partitions", {
  models <- c("GTR+G(4)+I", "HKY+G(4)", "GTR+G(4)+I")
  out <- evomodelCmds(model.vector = models)

  expect_named(out, c("grouped_commands", "all_commands", "partition_map", "raw_data"))
  # Partitions 1 and 3 share the same model, so should be grouped together
  expect_equal(nrow(out$grouped_commands), 2)
  expect_true("1,3" %in% out$grouped_commands$partitions)
  expect_true("2" %in% out$grouped_commands$partitions)
})

test_that("evomodelCmds produces valid lset/prset command syntax", {
  out <- evomodelCmds(model.vector = c("GTR+G", "HKY"))

  expect_true(any(grepl("^lset applyto=\\(1\\) nst=6 rates=gamma", out$all_commands)))
  expect_true(any(grepl("^lset applyto=\\(2\\) nst=2 rates=equal", out$all_commands)))
  expect_true(any(grepl("^prset applyto=\\(1\\)", out$all_commands)))
})

test_that("evomodelCmds partition_map reflects the original model per position", {
  out <- evomodelCmds(model.vector = c("JC", "GTR+I"))

  expect_equal(unname(out$partition_map), c("JC", "GTR+I"))
  expect_equal(names(out$partition_map), c("1", "2"))
})

test_that("evomodelCmds defaults an unrecognized model to GTR with a warning", {
  expect_warning(out <- evomodelCmds(model.vector = "NOTAMODEL"), "not recognized")
  expect_equal(out$raw_data$nst, 6)
})

test_that("evomodelCmds errors on empty or non-character input", {
  expect_error(evomodelCmds(model.vector = character(0)), "non-empty character")
  expect_error(evomodelCmds(model.vector = 123), "non-empty character")
})

test_that("evomodelCmds handles invariant-sites-only models", {
  out <- evomodelCmds(model.vector = "HKY+I")
  expect_equal(out$raw_data$rates, "propinv")
})

test_that("evomodelCmds keeps SYM and GTR apart through their base frequencies", {
  out <- evomodelCmds(model.vector = c("GTR+G(4)+I", "HKY+G(4)", "SYM+G(4)+I"))

  expect_equal(nrow(out$grouped_commands), 3)
  expect_true("lset applyto=(3) nst=6 rates=invgamma ngammacat=4;" %in% out$all_commands)
  expect_true("prset applyto=(1) statefreqpr=dirichlet(1,1,1,1);" %in% out$all_commands)
  expect_true("prset applyto=(3) statefreqpr=fixed(equal);" %in% out$all_commands)
  expect_equal(out$grouped_commands$model, c("GTR+G(4)+I", "HKY+G(4)", "SYM+G(4)+I"))
})

test_that("evomodelCmds never fixes the substitution rates", {
  out <- evomodelCmds(model.vector = c("TVM+G(4)", "TIM2e+I", "K80", "GTR"))
  expect_false(any(grepl("revmatpr", out$all_commands)))
})

test_that("evomodelCmds maps models to the right MrBayes scheme and frequencies", {
  out <- evomodelCmds(model.vector = c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TVMe", "TVM"))
  raw <- out$raw_data

  expect_equal(raw$nst, c(1, 1, 2, 2, 6, 6, 6, 6))
  expect_equal(raw$statefreqpr == "fixed(equal)",
               c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
})

test_that("evomodelCmds flags models MrBayes can only approximate", {
  out <- evomodelCmds(model.vector = c("GTR+G(4)", "TVM+G(4)"))

  expect_true(any(grepl("^\\[ Partition 2 \\(TVM\\+G\\(4\\)\\): not available in MrBayes; nst=6",
                        out$all_commands)))
  expect_false(any(grepl("Partition 1", out$all_commands)))
  # TVM is run as GTR in MrBayes, so both share one group named after both models
  expect_equal(out$grouped_commands$partitions, "1,2")
  expect_equal(out$grouped_commands$model, "GTR+G(4) / TVM+G(4)")
})
