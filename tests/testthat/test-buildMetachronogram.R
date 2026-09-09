# .graft_clade()'s correctness (exact target age, untouched-part invariance)
# is validated end-to-end against a real published dataset (Ringelberg et al.,
# reproducing the original Python/Phyx "graft_crown_clades.py" pipeline
# bit-for-bit on tip set and to floating-point precision on branch lengths);
# these unit tests cover the catGenes-facing interface and its edge cases
# with small synthetic trees.

.metachrono_backbone <- function() {
  set.seed(3)
  tr <- ape::rcoal(8, tip.label = paste0("B", 1:8))
  tr$edge.length <- tr$edge.length * 10
  tr
}

.metachrono_subtree <- function(n = 5, scale = 2, seed = 11, prefix = "S") {
  set.seed(seed)
  tr <- ape::rcoal(n, tip.label = paste0(prefix, 1:n))
  tr$edge.length <- tr$edge.length * scale
  tr
}

test_that("buildMetachronogram grafts a clade at exactly its target age", {
  backbone <- .metachrono_backbone()
  subtree <- .metachrono_subtree()

  bt <- ape::branching.times(backbone)
  cand_node <- NULL
  for (nd in as.integer(names(bt))) {
    if (ape::Ntip(ape::extract.clade(backbone, nd)) == 3) {
      cand_node <- nd
      break
    }
  }
  mrca_tips <- ape::extract.clade(backbone, cand_node)$tip.label
  target_age <- unname(bt[as.character(cand_node)])

  result <- buildMetachronogram(
    backbone = backbone,
    grafts = list(list(mrca = mrca_tips, subtree = subtree)),
    save = FALSE, verbose = FALSE
  )

  expect_s3_class(result, "phylo")
  expect_equal(ape::Ntip(result), ape::Ntip(backbone) - length(mrca_tips) + ape::Ntip(subtree))
  expect_true(all(subtree$tip.label %in% result$tip.label))
  # the replaced clade's own (backbone) tip names must be gone (subtree tips
  # use a disjoint "S" prefix, so no coincidental name reuse is possible here)
  expect_false(any(mrca_tips %in% result$tip.label))

  new_age <- unname(ape::branching.times(result)[as.character(ape::getMRCA(result, subtree$tip.label))])
  expect_equal(new_age, target_age, tolerance = 1e-8)
})

test_that("buildMetachronogram leaves the untouched part of the tree unchanged", {
  backbone <- .metachrono_backbone()
  subtree <- .metachrono_subtree()
  mrca_tips <- c("B1", "B2", "B4")

  result <- buildMetachronogram(
    backbone = backbone,
    grafts = list(list(mrca = mrca_tips, subtree = subtree)),
    save = FALSE, verbose = FALSE
  )

  kept_tips <- setdiff(backbone$tip.label, mrca_tips)
  d_before <- ape::cophenetic.phylo(ape::keep.tip(backbone, kept_tips))
  d_after <- ape::cophenetic.phylo(ape::keep.tip(result, kept_tips))
  d_after <- d_after[rownames(d_before), colnames(d_before)]

  expect_equal(d_after, d_before, tolerance = 1e-8)
})

test_that("buildMetachronogram supports posterior samples (multiPhylo) with matched pairing", {
  backbone <- .metachrono_backbone()
  backbone_list <- rep(list(backbone), 3)
  class(backbone_list) <- "multiPhylo"

  subtrees <- lapply(1:3, function(i) .metachrono_subtree(seed = 10 + i))
  class(subtrees) <- "multiPhylo"

  result <- buildMetachronogram(
    backbone = backbone_list,
    grafts = list(list(mrca = c("B1", "B2", "B4"), subtree = subtrees)),
    save = FALSE, verbose = FALSE
  )

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 3)
  expect_true(all(sapply(result, ape::Ntip) == ape::Ntip(result[[1]])))
})

test_that("buildMetachronogram recycles a single subtree across a multiPhylo backbone", {
  backbone <- .metachrono_backbone()
  backbone_list <- rep(list(backbone), 3)
  class(backbone_list) <- "multiPhylo"
  subtree <- .metachrono_subtree()

  result <- buildMetachronogram(
    backbone = backbone_list,
    grafts = list(list(mrca = c("B1", "B2", "B4"), subtree = subtree)),
    save = FALSE, verbose = FALSE
  )

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 3)
})

test_that("buildMetachronogram supports hierarchical (nested) grafting in list order", {
  backbone <- .metachrono_backbone()
  outer_subtree <- .metachrono_subtree(n = 4, scale = 3, seed = 21, prefix = "S")
  inner_subtree <- .metachrono_subtree(n = 3, scale = 1, seed = 22, prefix = "I")

  result <- buildMetachronogram(
    backbone = backbone,
    grafts = list(
      list(mrca = c("B1", "B2", "B4"), subtree = outer_subtree),
      list(mrca = outer_subtree$tip.label[1:2], subtree = inner_subtree)
    ),
    save = FALSE, verbose = FALSE
  )

  expect_s3_class(result, "phylo")
  expect_true(all(inner_subtree$tip.label %in% result$tip.label))
  expect_true(ape::is.ultrametric(result, tol = 1e-6))
})

test_that("buildMetachronogram saves a Newick file when save = TRUE", {
  backbone <- .metachrono_backbone()
  subtree <- .metachrono_subtree()
  outdir <- withr::local_tempdir()

  buildMetachronogram(
    backbone = backbone,
    grafts = list(list(mrca = c("B1", "B2", "B4"), subtree = subtree)),
    save = TRUE, dir = outdir, filename = "test_metachrono", verbose = FALSE
  )

  foldername <- file.path(outdir, format(Sys.time(), "%d%b%Y"))
  expect_true(file.exists(file.path(foldername, "test_metachrono.tre")))
})

test_that("buildMetachronogram errors when backbone or grafts is missing", {
  backbone <- .metachrono_backbone()
  subtree <- .metachrono_subtree()

  expect_error(buildMetachronogram(backbone = NULL, grafts = list()), "backbone")
  expect_error(
    buildMetachronogram(backbone = backbone, grafts = NULL),
    "grafts"
  )
})

test_that("buildMetachronogram errors when a graft's mrca has fewer than two tips", {
  backbone <- .metachrono_backbone()
  subtree <- .metachrono_subtree()

  expect_error(
    buildMetachronogram(
      backbone = backbone,
      grafts = list(list(mrca = "B1", subtree = subtree)),
      save = FALSE, verbose = FALSE
    ),
    "at least two"
  )
})

test_that("buildMetachronogram errors on a multiPhylo subtree of mismatched length", {
  backbone <- .metachrono_backbone()
  backbone_list <- rep(list(backbone), 3)
  class(backbone_list) <- "multiPhylo"

  subtrees <- lapply(1:2, function(i) .metachrono_subtree(seed = 30 + i))
  class(subtrees) <- "multiPhylo"

  expect_error(
    buildMetachronogram(
      backbone = backbone_list,
      grafts = list(list(mrca = c("B1", "B2", "B4"), subtree = subtrees)),
      save = FALSE, verbose = FALSE
    ),
    "backbone"
  )
})

test_that("buildMetachronogram errors when the replacement subtree has fewer than two tips", {
  backbone <- .metachrono_backbone()
  one_tip <- ape::drop.tip(.metachrono_subtree(n = 2), "S2")

  expect_error(
    buildMetachronogram(
      backbone = backbone,
      grafts = list(list(mrca = c("B1", "B2", "B4"), subtree = one_tip)),
      save = FALSE, verbose = FALSE
    ),
    "at least two tips"
  )
})
