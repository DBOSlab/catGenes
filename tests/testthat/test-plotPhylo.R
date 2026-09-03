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

test_that("a single highlight.color is recycled across every highlight.taxa", {
  # .col.tiplabels() expects tiplabels with underscores already converted to
  # spaces, matching the conversion plotPhylo() itself applies (R/plotPhylo.R,
  # `tree@phylo$tip.label <- gsub("_", " ", ...)`) before calling this helper.
  tiplabels <- c("Genus alpha", "Genus beta", "Genus gamma")

  tipdata <- catGenes:::.col.tiplabels(
    tiplabels = tiplabels,
    highlight.taxa = c("Genus_alpha", "Genus_beta"),
    highlight.color = "red"
  )

  expect_equal(tipdata$tocolor, c("red", "red", "black"))
})

test_that(".col.tiplabels still supports one distinct color per taxon", {
  tiplabels <- c("Genus alpha", "Genus beta", "Genus gamma")

  tipdata <- catGenes:::.col.tiplabels(
    tiplabels = tiplabels,
    highlight.taxa = c("Genus_alpha", "Genus_beta"),
    highlight.color = c("red", "blue")
  )

  expect_equal(tipdata$tocolor, c("red", "blue", "black"))
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

test_that("plotPhylo renders the full documented example (clades, gene labels, side phylogram)", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  data(Harpalyce_parsimony_tree, package = "catGenes", envir = environment())
  data(Harpalyce_raxml_tree, package = "catGenes", envir = environment())
  outdir <- withr::local_tempdir()

  Harpalyce_clade <- c("Harpalyce_brasiliana_Cardoso2510", "Harpalyce_formosa_Hughes2109")
  outgroup_taxa <- c("Dermatophyllum_secundiflorum", "Clathrotropis_nitida",
                     "Bowdichia_virgilioides")

  p <- suppressWarnings(
    plotPhylo(
      tree = Harpalyce_bayes_tree,
      layout = "rectangular",
      branch.width = 0.5,
      branch.supports = TRUE,
      add.raxml.tree = Harpalyce_raxml_tree,
      add.parsi.tree = Harpalyce_parsimony_tree,
      highlight.clade = Harpalyce_clade,
      fill.gradient = "#D53E4F",
      show.tip.label = TRUE,
      size.tip.label = 4,
      fontface.tip.label = "italic",
      understate.taxa = outgroup_taxa,
      gene.label = c("ITS/5.8S", "ETS", "matK", "trnL intron"),
      size.gene.label = 12,
      ylim.gene.label = NULL,
      phylogram.side = TRUE,
      phylogram.supports = TRUE,
      phylogram.height = 25,
      save = TRUE,
      dir = outdir,
      format = "pdf"
    )
  )

  expect_s3_class(p, "ggtree")

  foldername <- file.path(
    outdir,
    format(Sys.time(), "%d%b%Y")
  )

  expect_true(
    file.exists(
      file.path(foldername, "edited_tree_rectangular.pdf")
    ))
})

test_that("plotPhylo replaces taxon names and disables the support legend", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  tips <- Harpalyce_bayes_tree@phylo$tip.label

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 replace.taxa = setNames("Renamed_taxon", tips[1]),
                 support.legend = FALSE,
                 branch.supports = FALSE,
                 save = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plotPhylo custom xlim.tree/hexpand render without error", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 xlim.tree = c(0, 1),
                 hexpand = 0.2,
                 save = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plotPhylo fancy.tip.label renders Markdown-formatted tip labels with highlighted taxa", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  tips <- Harpalyce_bayes_tree@phylo$tip.label

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 fancy.tip.label = TRUE,
                 highlight.taxa = tips[1:2],
                 highlight.color = "blue",
                 save = FALSE)

  expect_s3_class(p, "ggtree")
})

test_that("plotPhylo fancy.tip.label renders with understated taxa", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  tips <- Harpalyce_bayes_tree@phylo$tip.label

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 fancy.tip.label = TRUE,
                 understate.taxa = tips[3:4],
                 save = FALSE)

  expect_s3_class(p, "ggtree")
})

test_that("plotPhylo fancy.tip.label combines highlight + understate + abbreviated genus", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  tips <- Harpalyce_bayes_tree@phylo$tip.label

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 fancy.tip.label = TRUE,
                 highlight.taxa = tips[1:2],
                 highlight.color = "blue",
                 understate.taxa = tips[3:4],
                 abbrev.tip.label = TRUE,
                 save = FALSE)

  expect_s3_class(p, "ggtree")
})

test_that("plotPhylo saves to PNG, JPG, and TIFF formats", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  outdir <- withr::local_tempdir()

  for (fmt in c("png", "jpg", "tiff")) {
    plotPhylo(tree = Harpalyce_bayes_tree, save = TRUE, dir = outdir,
              filename = paste0("tree_", fmt), format = fmt)
  }

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "tree_png.png")))
  expect_true(file.exists(file.path(foldername, "tree_jpg.jpg")))
  expect_true(file.exists(file.path(foldername, "tree_tiff.tiff")))
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
