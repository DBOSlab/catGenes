#' Build a metachronogram by grafting dated subtrees onto a backbone chronogram
#'
#' @author Domingos Cardoso
#'
#' @description Builds a large, composite time-calibrated phylogeny (a
#' "metachronogram", sensu Ringelberg et al.) by grafting one or more densely
#' sampled, ultrametric subtrees onto specified clades of a backbone
#' chronogram, replacing each targeted clade with the corresponding subtree
#' rescaled to exactly the age of the clade it replaces. This lets you
#' combine phylogenetic datasets generated at different taxonomic scales or
#' from different sequencing efforts (e.g. grafting densely sampled,
#' Sanger-sequenced genus-level trees onto a sparser, genome-scale backbone)
#' into a single time-scaled tree, while keeping every untouched part of the
#' backbone exactly as it was.
#'
#' @usage
#' buildMetachronogram(backbone = NULL,
#'                     grafts = NULL,
#'                     verbose = TRUE,
#'                     save = TRUE,
#'                     dir = "RESULTS_buildMetachronogram",
#'                     filename = "metachronogram")
#'
#' @param backbone An ultrametric \code{phylo} tree, or a \code{multiPhylo}
#' object (e.g. a posterior sample of backbone trees), onto which the
#' subtrees in \code{grafts} will be grafted.
#'
#' @param grafts A list of grafts to apply, in order. Each element is itself
#' a list with two components: \code{mrca}, a vector of two or more tip
#' labels whose most recent common ancestor in the (current state of the)
#' backbone defines the clade to be replaced; and \code{subtree}, the
#' ultrametric replacement \code{phylo} tree (or a \code{multiPhylo} of the
#' same length as \code{backbone}, for grafting a matching posterior subtree
#' sample onto each backbone tree). Grafts are applied in the order given, so
#' for hierarchical grafting (grafting a subtree inside a clade that was
#' itself just grafted in), list the least nested graft first.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each graft as
#' it is applied will not be printed in the console.
#'
#' @param save Logical, if \code{FALSE} the resulting metachronogram(s) are
#' returned as an R object but not saved to disk.
#'
#' @param dir The path to the directory where the metachronogram should be
#' saved. The default is to create a directory named
#' **RESULTS_buildMetachronogram** and save the tree within a subfolder named
#' after the current date.
#'
#' @param filename Name (without extension) of the saved Newick tree file.
#'
#' @details
#' For each graft, the clade defined by \code{mrca} is located in the
#' backbone and its age (time before the present, assuming the backbone tips
#' are all at time 0) is taken as the target age. The replacement subtree is
#' rescaled proportionally so that its own root height exactly equals that
#' target age, and is then grafted in at precisely the point the original
#' clade occupied, preserving the stem branch leading to it. The rest of the
#' backbone — including branch lengths and node ages outside the replaced
#' clade — is left unchanged.
#'
#' @return A \code{phylo} object (if \code{backbone} is a single tree) or a
#' \code{multiPhylo} object (if \code{backbone} is a \code{multiPhylo}), with
#' every requested graft applied.
#'
#' @seealso \code{\link{plotChronogram}}
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' # single backbone tree, single graft
#' metatree <- buildMetachronogram(
#'   backbone = my_backbone_tree,
#'   grafts = list(
#'     list(mrca = c("Genus_speciesA", "Genus_speciesB"),
#'         subtree = my_genus_level_tree)
#'   )
#' )
#'
#' # posterior sample of backbone trees, two grafts applied to every tree,
#' # each paired with its own posterior sample of the replacement subtree
#' metatrees <- buildMetachronogram(
#'   backbone = my_backbone_posterior,
#'   grafts = list(
#'     list(mrca = c("Genus1_speciesA", "Genus1_speciesB"),
#'         subtree = genus1_posterior),
#'     list(mrca = c("Genus2_speciesA", "Genus2_speciesB"),
#'         subtree = genus2_posterior)
#'   ),
#'   dir = "RESULTS_buildMetachronogram"
#' )
#'}
#'
#' @importFrom ape getMRCA branching.times extract.clade node.depth.edgelength
#' @importFrom ape drop.tip bind.tree write.tree Ntip
#'
#' @export
#'
buildMetachronogram <- function(backbone = NULL,
                                grafts = NULL,
                                verbose = TRUE,
                                save = TRUE,
                                dir = "RESULTS_buildMetachronogram",
                                filename = "metachronogram") {

  if (is.null(backbone)) {
    stop("Please provide 'backbone', an ultrametric tree (or multiPhylo of ",
        "trees) to graft onto.")
  }
  if (is.null(grafts) || length(grafts) == 0) {
    stop("Please provide 'grafts', a list of MRCA + replacement subtree ",
        "specifications.")
  }
  if (!inherits(backbone, "phylo") && !inherits(backbone, "multiPhylo")) {
    stop("'backbone' must be a 'phylo' or 'multiPhylo' object.")
  }

  single_backbone <- inherits(backbone, "phylo")
  backbone_list <- if (single_backbone) list(backbone) else backbone
  n <- length(backbone_list)

  for (i in seq_along(grafts)) {
    g <- grafts[[i]]
    if (is.null(g$mrca) || length(g$mrca) < 2) {
      stop("Graft ", i, ": 'mrca' must be a vector of at least two tip labels.")
    }
    if (is.null(g$subtree)) {
      stop("Graft ", i, ": 'subtree' must be provided.")
    }
    if (inherits(g$subtree, "phylo")) {
      grafts[[i]]$subtree <- rep(list(g$subtree), n)
    } else if (inherits(g$subtree, "multiPhylo")) {
      if (length(g$subtree) != n) {
        stop("Graft ", i, ": 'subtree' has ", length(g$subtree),
            " tree(s) but 'backbone' has ", n, "; provide either a single ",
            "subtree (reused for every backbone tree) or exactly ", n, ".")
      }
      grafts[[i]]$subtree <- as.list(g$subtree)
    } else {
      stop("Graft ", i, ": 'subtree' must be a 'phylo' or 'multiPhylo' object.")
    }
  }

  results <- vector("list", n)

  for (i in seq_len(n)) {
    if (verbose) {
      message("Building metachronogram ", i, " of ", n, "...")
    }
    tree <- backbone_list[[i]]
    for (j in seq_along(grafts)) {
      if (verbose) {
        message("  Grafting clade ", j, " of ", length(grafts),
               " (MRCA: ", paste(grafts[[j]]$mrca, collapse = ", "), ")")
      }
      tree <- .graft_clade(tree, grafts[[j]]$mrca, grafts[[j]]$subtree[[i]])
    }
    results[[i]] <- tree
  }

  if (save) {
    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername, recursive = TRUE)
    }
    outfile <- paste0(foldername, "/", filename, ".tre")

    to_write <- results
    class(to_write) <- "multiPhylo"
    ape::write.tree(to_write, file = outfile)

    if (verbose) {
      message("\nSaved metachronogram", if (n > 1) "s" else "", " to: ", outfile)
    }
  }

  class(results) <- "multiPhylo"
  if (single_backbone) {
    return(invisible(results[[1]]))
  }
  return(invisible(results))
}


#-------------------------------------------------------------------------------
# Auxiliary function to graft a single rescaled subtree into a backbone tree,
# replacing the clade defined by mrca_tips with subtree, rescaled so its root
# height exactly equals the age of the clade it replaces. See buildMetachronogram().

.graft_clade <- function(tree, mrca_tips, subtree) {

  mrca_tips <- intersect(mrca_tips, tree$tip.label)
  if (length(mrca_tips) < 2) {
    stop("At least two of the tip labels defining an MRCA must be present ",
        "in the tree at the point this graft is applied.")
  }

  mrca_node <- ape::getMRCA(tree, mrca_tips)

  target_age <- unname(ape::branching.times(tree)[as.character(mrca_node)])

  old_clade <- ape::extract.clade(tree, mrca_node)

  if (ape::Ntip(subtree) < 2) {
    stop("The replacement subtree must have at least two tips.")
  }

  tip_depths <- ape::node.depth.edgelength(subtree)[seq_len(ape::Ntip(subtree))]
  root_height <- max(tip_depths)
  if (!is.finite(root_height) || root_height <= 0) {
    stop("The replacement subtree has no resolvable root height; check its ",
        "branch lengths.")
  }
  subtree$edge.length <- subtree$edge.length * (target_age / root_height)

  # Drop all but one tip of the old clade to leave a single placeholder tip
  # at the exact attachment point, graft the rescaled subtree onto it at the
  # target age, then remove the now-redundant placeholder.
  placeholder <- old_clade$tip.label[1]
  tips_to_drop <- setdiff(old_clade$tip.label, placeholder)
  pruned <- ape::drop.tip(tree, tips_to_drop)

  grafted <- ape::bind.tree(
    pruned, subtree,
    where = which(pruned$tip.label == placeholder),
    position = target_age
  )
  grafted <- ape::drop.tip(grafted, placeholder)

  return(grafted)
}
