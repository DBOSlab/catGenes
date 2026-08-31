# plotPhylo() renders fast (< 1s) even on the full ~87-tip Harpalyce tree, so
# these tests exercise the real function against the package's bundled
# treedata objects rather than mocking ggtree/phytools internals.

test_that("plotPhylo returns a ggtree/ggplot object for a single tree", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())

  p <- plotPhylo(tree = Harpalyce_bayes_tree, save = FALSE)

  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "ggtree")
})

test_that("plotPhylo accepts alternative layouts", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())

  for (lay in c("rectangular", "circular", "fan")) {
    p <- plotPhylo(tree = Harpalyce_bayes_tree, layout = lay, save = FALSE)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plotPhylo merges node support from additional RAxML/parsimony trees", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  data(Harpalyce_raxml_tree, package = "catGenes", envir = environment())
  data(Harpalyce_parsimony_tree, package = "catGenes", envir = environment())

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 add.raxml.tree = Harpalyce_raxml_tree,
                 add.parsi.tree = Harpalyce_parsimony_tree,
                 save = FALSE)

  expect_s3_class(p, "ggtree")
})

test_that("plotPhylo highlights specified taxa without erroring", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  tips <- Harpalyce_bayes_tree@phylo$tip.label

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 highlight.taxa = tips[1:3],
                 highlight.color = "red",
                 save = FALSE)

  expect_s3_class(p, "ggtree")
})

test_that("plotPhylo prunes the requested taxa out of the tree", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  tips <- Harpalyce_bayes_tree@phylo$tip.label

  p <- suppressMessages(plotPhylo(tree = Harpalyce_bayes_tree,
                                  prune.taxa = tips[1:3],
                                  save = FALSE))

  expect_s3_class(p, "ggplot")
  remaining_labels <- p$data$label[p$data$isTip]
  expect_false(any(tips[1:3] %in% remaining_labels))
})

test_that("plotPhylo abbreviates tip labels via abbrevGen when requested", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())

  p <- plotPhylo(tree = Harpalyce_bayes_tree, abbrev.tip.label = TRUE, save = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plotPhylo saves a file to disk when save = TRUE", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  outdir <- withr::local_tempdir()

  plotPhylo(tree = Harpalyce_bayes_tree,
           save = TRUE,
           dir = outdir,
           filename = "mytree",
           format = "pdf")

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "mytree.pdf")))
})
