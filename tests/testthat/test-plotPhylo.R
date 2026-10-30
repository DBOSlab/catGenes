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

test_that("collapse.clade turns named clades into triangles labelled with the clade name", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  mex <- c("Harpalyce_arborescens_Lezama19187", "Harpalyce_sp6_Nunez550")
  poe <- c("Poecilanthe_sp2_Saddi933", "Poecilanthe_ulei_Queiroz9584")

  p <- plotPhylo(tree = Harpalyce_bayes_tree,
                 collapse.clade = list("Mexican clade" = mex, poe),
                 save = FALSE)

  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geoms == "GeomPolygon"), 2)

  texts <- unlist(lapply(p$layers[geoms == "GeomText"],
                         function(l) l$aes_params$label))
  expect_true("Mexican clade" %in% texts)

  # Terminals of the named clade are hidden; those of the unnamed one are kept
  d <- p$data[p$data$isTip, ]
  expect_true(all(is.na(d$x[grepl("^Harpalyce (arborescens|sp6|sp1) ", d$label)])))
  expect_false(any(is.na(d$x[grepl("^Poecilanthe ", d$label)])))
})

test_that("collapse.scale shrinks named triangles and keeps the plot valid", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  mex <- list("Mexican clade" = c("Harpalyce_arborescens_Lezama19187",
                                  "Harpalyce_sp6_Nunez550"))

  full <- plotPhylo(tree = Harpalyce_bayes_tree, collapse.clade = mex, save = FALSE)
  small <- plotPhylo(tree = Harpalyce_bayes_tree, collapse.clade = mex,
                     collapse.scale = 0.25, save = FALSE)

  triangle_height <- function(p) {
    poly <- Filter(function(l) inherits(l$geom, "GeomPolygon"), p$layers)[[1]]
    diff(range(poly$data$y))
  }
  expect_equal(triangle_height(small), triangle_height(full) * 0.25,
               tolerance = 1e-6)
})

test_that("collapse.clade works with vertical names and circular layouts", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  mex <- list("Mexican clade" = c("Harpalyce_arborescens_Lezama19187",
                                  "Harpalyce_sp6_Nunez550"))

  expect_s3_class(plotPhylo(tree = Harpalyce_bayes_tree, collapse.clade = mex,
                            collapse.label.orientation = "vertical",
                            save = FALSE), "ggtree")
  expect_s3_class(plotPhylo(tree = Harpalyce_bayes_tree, layout = "circular",
                            collapse.clade = mex, save = FALSE), "ggtree")
})

test_that("collapse.clade validates its input", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  mex <- c("Harpalyce_arborescens_Lezama19187", "Harpalyce_sp6_Nunez550")

  expect_error(plotPhylo(tree = Harpalyce_bayes_tree,
                         collapse.clade = c("Harpalyce_arborescens_Lezama19187", "Nope_nope"),
                         save = FALSE), "not found in the tree: Nope nope")
  expect_error(plotPhylo(tree = Harpalyce_bayes_tree,
                         collapse.clade = list(mex[1]), save = FALSE),
               "exactly two tip labels")
  expect_error(plotPhylo(tree = Harpalyce_bayes_tree,
                         collapse.clade = list(mex, c("Harpalyce_sp6_Nunez550",
                                                      "Harpalyce_goldmanii")),
                         save = FALSE), "Nested clades")
  expect_error(plotPhylo(tree = Harpalyce_bayes_tree, collapse.clade = mex,
                         collapse.label.orientation = "diagonal", save = FALSE),
               "horizontal' or 'vertical")
  expect_error(plotPhylo(tree = Harpalyce_bayes_tree, collapse.clade = mex,
                         collapse.scale = 0, save = FALSE), "collapse.scale")
})

test_that("size.tip.label is applied when tips are highlighted or understated", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())

  tip_sizes <- function(p) {
    unlist(lapply(p$layers, function(l) {
      if (inherits(l$geom, "GeomText")) l$aes_params$size
    }))
  }

  for (args in list(list(highlight.taxa = "Harpalyce", highlight.color = "red"),
                    list(understate.taxa = "Leptolobium"))) {
    p <- do.call(plotPhylo, c(list(tree = Harpalyce_bayes_tree,
                                   size.tip.label = 6,
                                   branch.supports = FALSE,
                                   save = FALSE), args))
    expect_true(6 %in% tip_sizes(p))
    expect_false(3 %in% tip_sizes(p))
  }
})

test_that("collapsing clades keeps the saved figure size and side phylogram placement", {
  data(Harpalyce_bayes_tree, package = "catGenes", envir = environment())
  mex <- list("Mexican clade" = c("Harpalyce_arborescens_Lezama19187",
                                  "Harpalyce_sp6_Nunez550"))
  saved <- list()
  testthat::local_mocked_bindings(
    .save_tree = function(fulltree, fullname, dpi, format, base_height, base_width) {
      saved[[length(saved) + 1]] <<- list(fulltree = fulltree,
                                          size = c(base_height, base_width))
    }
  )
  dir <- withr::local_tempdir()

  plotPhylo(tree = Harpalyce_bayes_tree, phylogram.side = TRUE,
            save = TRUE, dir = dir)
  plotPhylo(tree = Harpalyce_bayes_tree, phylogram.side = TRUE,
            collapse.clade = mex, save = TRUE, dir = dir)

  expect_false(anyNA(saved[[2]]$size))
  expect_equal(saved[[2]]$size, saved[[1]]$size)

  inset <- function(p) {
    l <- Filter(function(l) inherits(l$geom, "GeomCustomAnn"), p$layers)[[1]]
    unlist(l$geom_params[c("xmin", "xmax", "ymin", "ymax")])
  }
  expect_false(anyNA(inset(saved[[2]]$fulltree)))
  expect_equal(inset(saved[[2]]$fulltree), inset(saved[[1]]$fulltree))
})

test_that("decimal supports with an unsupported root stay on the 0-1 scale", {
  # RAxML-NG support trees (e.g. TBE) store supports as proportions and leave
  # the root without a label
  tr <- ape::read.tree(text = "((A:1,B:1)0.962714:1,((C:1,D:1)0.706500:1,E:1)0.891222:1);")
  tr <- treeio::as.treedata(tr)

  p <- plotPhylo(tree = tr, save = FALSE)
  d <- p$data[!p$data$isTip, ]

  expect_equal(sort(d$prob_numeric), c(0.7065, 0.891222, 0.962714, 1))
  expect_setequal(as.character(d$support_cat),
                  c(">= 1.0", "0.96 - 0.99", "0.86 - 0.95", "0.70 - 0.85"))
})

test_that("bootstrap supports with an unsupported root stay on the 0-100 scale", {
  tr <- ape::read.tree(text = "((A:1,B:1)96:1,((C:1,D:1)72:1,E:1)88:1);")
  tr <- treeio::as.treedata(tr)

  p <- plotPhylo(tree = tr, save = FALSE)
  d <- p$data[!p$data$isTip, ]

  expect_equal(sort(d$prob_numeric), c(72, 88, 96, 100))
  expect_setequal(as.character(d$support_cat), c("100", "96 - 99", "86 - 95", "70 - 85"))
})

test_that(".support_numeric reads combined and non-numeric supports", {
  expect_equal(catGenes:::.support_numeric(c("0.95/85/90", "100", "-", "", NA)),
               c(0.95, 100, NA, NA, NA))
})

test_that("support categories match their legend labels at the boundaries", {
  cat_of <- function(newick) {
    p <- plotPhylo(tree = treeio::as.treedata(ape::read.tree(text = newick)),
                   save = FALSE)
    d <- p$data[!p$data$isTip & !is.na(p$data$support_cat), ]
    setNames(as.character(d$support_cat), as.character(d$prob_numeric))
  }

  pp <- cat_of("(((((A:1,B:1)0.995:1,C:1)0.95:1,D:1)0.85:1,E:1)0.70:1,(F:1,G:1)0.69:1);")
  expect_equal(unname(pp[c("0.995", "0.95", "0.85", "0.7", "0.69")]),
               c("0.96 - 0.99", "0.86 - 0.95", "0.70 - 0.85", "0.70 - 0.85", "< 0.70"))

  bs <- cat_of("(((((A:1,B:1)99.6:1,C:1)95:1,D:1)85:1,E:1)70:1,(F:1,G:1)69:1);")
  expect_equal(unname(bs[c("99.6", "95", "85", "70", "69")]),
               c("96 - 99", "86 - 95", "70 - 85", "70 - 85", "< 70"))
})

test_that("decimal support labels show two decimals, truncated like their category", {
  tr <- treeio::as.treedata(ape::read.tree(text = "(((A:1,B:1)0.996:1,C:1)0.9:1,D:1);"))
  p <- plotPhylo(tree = tr, save = FALSE)

  text_layers <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomText"),
                              logical(1)))
  labs <- unlist(lapply(text_layers, function(i) ggplot2::layer_data(p, i)$label))
  expect_true(all(c("0.99", "0.90") %in% labs))
  expect_false(any(c("1.00", "0.9") %in% labs))
})
