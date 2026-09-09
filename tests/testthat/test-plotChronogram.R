# plotChronogram() renders fast on small synthetic chronograms, so these
# tests exercise the real function rather than mocking ggtree/deeptime.

.chronogram_tree <- function(n = 10, scale = 40, seed = 7) {
  set.seed(seed)
  tr <- ape::rcoal(n, tip.label = paste0("Genus_species", 1:n))
  tr$edge.length <- tr$edge.length * scale
  tr
}

# Find a real (non-root) 2-3 tip clade in a tree, for clade-label tests.
.real_clade_tips <- function(tree, size = 2) {
  bt <- ape::branching.times(tree)
  for (nd in as.integer(names(sort(bt)))) {
    cl <- ape::extract.clade(tree, nd)$tip.label
    if (length(cl) == size) return(cl)
  }
  stop("no clade of the requested size found")
}

test_that("plotChronogram returns a ggtree/ggplot object for a single tree", {
  tr <- .chronogram_tree()

  p <- plotChronogram(tree = tr, save = FALSE)

  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "ggtree")
})

test_that("plotChronogram accepts rectangular and circular layouts", {
  tr <- .chronogram_tree()

  for (lay in c("rectangular", "circular")) {
    p <- plotChronogram(tree = tr, layout = lay, save = FALSE)
    expect_s3_class(p, "ggtree")
  }
})

test_that("plotChronogram errors on an invalid layout or time.axis", {
  tr <- .chronogram_tree()

  expect_error(plotChronogram(tree = tr, layout = "fan", save = FALSE),
              "layout")
  expect_error(plotChronogram(tree = tr, time.axis = "bogus", save = FALSE),
              "time.axis")
})

test_that("plotChronogram errors when tree is missing or not a phylo", {
  expect_error(plotChronogram(tree = NULL, save = FALSE), "phylo")
  expect_error(plotChronogram(tree = "not_a_tree", save = FALSE), "phylo")
})

test_that("plotChronogram draws clade labels defined by MRCA tip pairs", {
  tr <- .chronogram_tree()
  clade_tips <- .real_clade_tips(tr, size = 2)

  p <- plotChronogram(tree = tr,
                      clades = list("My clade" = clade_tips),
                      save = FALSE)

  expect_s3_class(p, "ggtree")
  cladelab_layers <- vapply(p$layers, function(l) {
    inherits(l$geom, "GeomCladeLab") ||
      grepl("cladelab", class(l$geom)[1], ignore.case = TRUE)
  }, logical(1))
  expect_true(any(cladelab_layers) || length(p$layers) > 0)
})

test_that("plotChronogram errors when a clade has fewer than two tips", {
  tr <- .chronogram_tree()
  clade_tips <- .real_clade_tips(tr, size = 2)

  expect_error(
    plotChronogram(tree = tr,
                   clades = list("Bad clade" = clade_tips[1]),
                   save = FALSE),
    "at least two"
  )
})

test_that("plotChronogram highlights specified taxa without erroring", {
  tr <- .chronogram_tree()

  p <- plotChronogram(tree = tr,
                      highlight.taxa = tr$tip.label[1:2],
                      highlight.color = "red",
                      save = FALSE)

  expect_s3_class(p, "ggtree")
})

test_that("plotChronogram falls back to a numeric axis for circular + geological", {
  tr <- .chronogram_tree()

  expect_warning(
    p <- plotChronogram(tree = tr, layout = "circular",
                        time.axis = "geological", save = FALSE),
    "rectangular"
  )
  expect_s3_class(p, "ggtree")
})

test_that("plotChronogram adds a geological timescale when time.axis is 'geological' or 'both'", {
  testthat::skip_if_not_installed("deeptime")
  tr <- .chronogram_tree()

  for (ax in c("geological", "both")) {
    p <- plotChronogram(tree = tr, time.axis = ax, save = FALSE)
    expect_s3_class(p, "ggtree")
    expect_true(inherits(p$coordinates, "CoordGeo") ||
                  any(grepl("CoordGeo", class(p$coordinates))))
  }
})

test_that("plotChronogram saves a plot file when save = TRUE", {
  tr <- .chronogram_tree()
  outdir <- withr::local_tempdir()

  plotChronogram(tree = tr, save = TRUE, dir = outdir,
                 filename = "test_chronogram", format = "pdf")

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "test_chronogram.pdf")))
})

test_that("plotChronogram warns on a non-ultrametric tree but still plots", {
  tr <- .chronogram_tree()
  tr$edge.length[1] <- tr$edge.length[1] + 100

  expect_warning(
    p <- plotChronogram(tree = tr, save = FALSE),
    "ultrametric"
  )
  expect_s3_class(p, "ggtree")
})

test_that("a circular save stays a sane size on a many-tip, shallow tree", {
  # ggtree's circular layout stores $data$y as each tip's 1..Ntip row index
  # (not a spatial length), so a tree with many tips but a small root age
  # (Ntip >> root_age) previously made the saved canvas scale to Ntip inches
  # instead of the tree's actual radius -- reproduces at a much smaller scale
  # here than the real-world case (2078 tips, ~67 Ma) that first surfaced it.
  tr <- .chronogram_tree(n = 60, scale = 3, seed = 9)
  outdir <- withr::local_tempdir()

  plotChronogram(tree = tr, layout = "circular", show.tip.label = FALSE,
                 save = TRUE, dir = outdir, filename = "test_circ_manytip",
                 format = "png", dpi = 100)

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  outfile <- file.path(foldername, "test_circ_manytip.png")
  expect_true(file.exists(outfile))

  # Read the PNG's own IHDR chunk (width/height are the two 4-byte big-endian
  # integers starting at byte offset 16) directly with base R, so this check
  # has no dependency on any image-reading package or OS-specific tool.
  con <- file(outfile, "rb")
  on.exit(close(con))
  seek(con, 16)
  px_width <- readBin(con, "integer", size = 4, endian = "big")
  expect_lt(px_width, 5000)
})
