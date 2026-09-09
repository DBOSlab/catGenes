#' Plot a time-calibrated phylogeny (chronogram or metachronogram)
#'
#' @author Domingos Cardoso
#'
#' @description Function to plot a [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html)-based
#' time-calibrated phylogeny, such as a chronogram or a metachronogram built
#' with \code{\link{buildMetachronogram}}, using real branch lengths on a
#' calendar-time axis (present = 0, root = oldest). The x-axis can show plain
#' numeric time (in millions of years, unless the input tree is scaled
#' otherwise), a geological timescale (epochs, periods, eras and/or eons, via
#' the \href{https://cran.r-project.org/package=deeptime}{deeptime} package),
#' or both. Named clades can be labelled around the tree, each defined by the
#' most recent common ancestor (MRCA) of two or more of its tip labels.
#'
#' @usage
#' plotChronogram(tree = NULL,
#'                layout = "rectangular",
#'                time.axis = "numeric",
#'                geo.ranks = c("epochs", "periods"),
#'                clades = NULL,
#'                clade.label.offset = NULL,
#'                clade.label.offset.text = NULL,
#'                clade.label.fontsize = 3.5,
#'                clade.label.barsize = 1.2,
#'                highlight.taxa = NULL,
#'                highlight.color = NULL,
#'                show.tip.label = TRUE,
#'                size.tip.label = 3,
#'                fontface.tip.label = "italic",
#'                branch.width = 0.5,
#'                save = FALSE,
#'                dpi = 600,
#'                dir = "RESULTS_plotChronogram",
#'                filename = NULL,
#'                format = "pdf",
#'                ...)
#'
#' @param tree An ultrametric \code{phylo} object with branch lengths in
#' units of time (e.g. millions of years), such as the output of
#' \code{\link{buildMetachronogram}}.
#'
#' @param layout Either \code{"rectangular"} or \code{"circular"}.
#'
#' @param time.axis One of \code{"numeric"} (a plain numeric time axis,
#' millions of years before present), \code{"geological"} (a geological
#' timescale only, via \code{deeptime::coord_geo()}), or \code{"both"}
#' (numeric axis plus the geological timescale). \code{"geological"} and
#' \code{"both"} require the \pkg{deeptime} package and only apply when
#' \code{layout = "rectangular"}; a \code{"circular"} layout always uses a
#' plain numeric axis, regardless of \code{time.axis}.
#'
#' @param geo.ranks A character vector naming which geological ranks to show
#' as stacked bands when \code{time.axis} is \code{"geological"} or
#' \code{"both"}. Any combination of \code{"eons"}, \code{"eras"},
#' \code{"periods"}, \code{"epochs"} and \code{"stages"} (see
#' \code{deeptime::coord_geo()}). Defaults to \code{c("epochs", "periods")}.
#'
#' @param clades A named list defining clades to label around the tree. Each
#' name is the label to print, and each element is a vector of two or more
#' tip labels whose MRCA in \code{tree} defines that clade -- the same
#' convention used by the \code{mrca} element of \code{\link{buildMetachronogram}}'s
#' \code{grafts}. Defaults to \code{NULL}, i.e. no clade labels are drawn.
#'
#' @param clade.label.offset Horizontal offset of the clade label bracket
#' from the tips it labels. Defaults to \code{NULL}, which sets it
#' automatically to 15% of the tree's root age.
#'
#' @param clade.label.offset.text Horizontal offset of the clade label text
#' from its bracket. Defaults to \code{NULL}, which sets it automatically to
#' 2% of the tree's root age.
#'
#' @param clade.label.fontsize Font size of the clade labels. Defaults to 3.5.
#'
#' @param clade.label.barsize Line width of the clade label brackets.
#' Defaults to 1.2.
#'
#' @param highlight.taxa Define a vector with specific tip labels or the name of
#' a particular genus/species to be highlighted with a color in the tree.
#'
#' @param highlight.color Color vector to highlight specific tip labels in the
#' tree, as defined with \code{highlight.taxa}. A single color is recycled
#' across all elements of \code{highlight.taxa}.
#'
#' @param show.tip.label Logical, if \code{FALSE} the tip labels will not be
#' plotted.
#'
#' @param size.tip.label A numeric value for the size of tip labels. Defaults
#' to 3.
#'
#' @param fontface.tip.label The font face for tip labels: 1 (plain, default),
#' 2 (bold), 3 (italic), 4 (bold italic). Defaults to \code{"italic"}.
#'
#' @param branch.width A numeric vector giving the width of the branches of
#' the plotted phylogeny. Defaults to 0.5.
#'
#' @param save Logical, if \code{TRUE} the plotted chronogram will be saved
#' on disk.
#'
#' @param dpi A number specifying the plot resolution when \code{format} is a
#' raster format (\code{"png"}, \code{"tiff"} or \code{"jpg"}). Defaults to 600.
#'
#' @param dir The path to the directory where the plot should be saved. The
#' default is to create a directory named **RESULTS_plotChronogram** and save
#' the tree within a subfolder named after the current date.
#'
#' @param filename Name (without extension) of the saved plot. Defaults to
#' \code{"chronogram_<layout>"}.
#'
#' @param format Format of the saved plot: \code{"pdf"} (default), \code{"png"},
#' \code{"tiff"} or \code{"jpg"}.
#'
#' @param ... Other arguments passed on to \code{ggtree::geom_tiplab()}.
#'
#' @return A ggtree/ggplot object with the plotted chronogram, invisibly
#' saved to disk if \code{save = TRUE}.
#'
#' @seealso \code{\link{buildMetachronogram}}, \code{\link{plotPhylo}}
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' plotChronogram(tree = my_metachronogram,
#'               layout = "rectangular",
#'               time.axis = "both",
#'               clades = list("Clade A" = c("Genus_speciesA", "Genus_speciesB"),
#'                             "Clade B" = c("Genus_speciesC", "Genus_speciesD")),
#'               highlight.taxa = "Genus_speciesA",
#'               highlight.color = "red",
#'               save = TRUE,
#'               dir = "results_chronogram",
#'               filename = "my_chronogram",
#'               format = "pdf")
#'}
#'
#' @importFrom ggtree ggtree revts geom_tiplab geom_cladelab theme_tree2 "%<+%"
#' @importFrom ggplot2 aes scale_x_continuous theme margin
#' @importFrom cowplot save_plot
#' @importFrom ape getMRCA Ntip node.depth.edgelength is.ultrametric
#'
#' @export
#'
plotChronogram <- function(tree = NULL,
                           layout = "rectangular",
                           time.axis = "numeric",
                           geo.ranks = c("epochs", "periods"),
                           clades = NULL,
                           clade.label.offset = NULL,
                           clade.label.offset.text = NULL,
                           clade.label.fontsize = 3.5,
                           clade.label.barsize = 1.2,
                           highlight.taxa = NULL,
                           highlight.color = NULL,
                           show.tip.label = TRUE,
                           size.tip.label = 3,
                           fontface.tip.label = "italic",
                           branch.width = 0.5,
                           save = FALSE,
                           dpi = 600,
                           dir = "RESULTS_plotChronogram",
                           filename = NULL,
                           format = "pdf",
                           ...) {

  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop("Please provide 'tree', a 'phylo' object with branch lengths in ",
        "units of time.")
  }
  if (!layout %in% c("rectangular", "circular")) {
    stop("'layout' must be either 'rectangular' or 'circular'.")
  }
  if (!time.axis %in% c("numeric", "geological", "both")) {
    stop("'time.axis' must be one of 'numeric', 'geological' or 'both'.")
  }
  if (!ape::is.ultrametric(tree, tol = 1e-6)) {
    warning("'tree' does not appear to be ultrametric; the time axis may ",
           "not be meaningful.")
  }

  use_geo <- time.axis %in% c("geological", "both") && layout == "rectangular"
  if (time.axis %in% c("geological", "both") && layout == "circular") {
    warning("A geological timescale is only supported for layout = ",
           "'rectangular'; using a plain numeric axis instead.")
  }
  if (use_geo && !requireNamespace("deeptime", quietly = TRUE)) {
    warning("Package 'deeptime' is required for time.axis = '", time.axis,
           "'; using a plain numeric axis instead. Install it with ",
           "install.packages(\"deeptime\").")
    use_geo <- FALSE
  }

  #-----------------------------------------------------------------------------
  # Create a new directory to save the plot

  if (save) {
    if (is.null(filename)) {
      filename <- paste0("chronogram_", layout)
    }

    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    fullname <- paste0(foldername, "/", filename, ".", format)
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername)
    }
  }

  #-----------------------------------------------------------------------------
  # Plotting the figure

  root_age <- max(ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))])

  tree$tip.label <- gsub("_", " ", tree$tip.label)
  tiplabels <- tree$tip.label

  tree_plot <- ggtree::ggtree(tree, layout = layout, ladderize = TRUE,
                              linewidth = branch.width)
  tree_plot <- ggtree::revts(tree_plot)

  # Rough estimate (in inches) of tip label text width from font size and
  # the longest tip label; reused below both to size the default clade-label
  # offsets and to size the saved plot. root_age can range from a few Ma to
  # hundreds of Ma across catGenes use cases, while the space needed to clear
  # tip label text does not scale with tree age, so this stays age-independent.
  char_width <- size.tip.label * 0.022
  tip_label_extent <- root_age * 0.01 +
    if (show.tip.label) max(nchar(tiplabels), 1) * char_width else char_width

  # Add clade labels, defined by the MRCA of >= 2 tip labels (same convention
  # as the 'mrca' element of buildMetachronogram()'s grafts)
  if (!is.null(clades)) {
    if (is.null(clade.label.offset)) {
      clade.label.offset <- tip_label_extent + char_width * 2
    }
    if (is.null(clade.label.offset.text)) {
      clade.label.offset.text <- char_width * 2
    }
    for (clade_name in names(clades)) {
      clade_tips <- gsub("_", " ", clades[[clade_name]])
      if (length(clade_tips) < 2) {
        stop("Clade '", clade_name, "': at least two tip labels are ",
            "required to define its MRCA.")
      }
      clade_node <- ape::getMRCA(tree, clade_tips)
      tree_plot <- tree_plot +
        ggtree::geom_cladelab(node = clade_node, label = clade_name,
                              offset = clade.label.offset,
                              offset.text = clade.label.offset.text,
                              fontsize = clade.label.fontsize,
                              barsize = clade.label.barsize)
    }
  }

  # Add tip labels, optionally highlighting specific taxa
  if (show.tip.label) {
    if (is.null(highlight.taxa)) {
      tree_plot <- tree_plot +
        ggtree::geom_tiplab(size = size.tip.label,
                            fontface = fontface.tip.label,
                            offset = root_age * 0.01,
                            show.legend = FALSE, ...)
    } else {
      tipdata <- .col.tiplabels(tiplabels = tiplabels,
                                highlight.taxa = highlight.taxa,
                                highlight.color = highlight.color)
      cols <- c("black", highlight.color)

      if (layout == "circular") {
        tree_plot <- tree_plot %<+% tipdata +
          ggtree::geom_tiplab(offset = root_age * 0.01, size = size.tip.label,
                              fontface = fontface.tip.label,
                              aes(angle = angle, color = tocolor),
                              show.legend = FALSE, ...) +
          ggplot2::scale_colour_manual(values = cols, breaks = cols,
                                       guide = "none")
      } else {
        tree_plot <- tree_plot %<+% tipdata +
          ggtree::geom_tiplab(offset = root_age * 0.01, size = size.tip.label,
                              fontface = fontface.tip.label,
                              aes(color = tocolor),
                              show.legend = FALSE, ...) +
          ggplot2::scale_colour_manual(values = cols, breaks = cols,
                                       guide = "none")
      }
    }
  }

  # Time axis
  if (layout == "rectangular") {
    pad_right <- root_age * 0.3
    pad_left <- root_age * 0.02

    if (time.axis %in% c("numeric", "both")) {
      tree_plot <- tree_plot +
        ggtree::theme_tree2() +
        ggplot2::scale_x_continuous(labels = abs)
    }
    if (use_geo) {
      tree_plot <- tree_plot +
        deeptime::coord_geo(dat = as.list(geo.ranks), neg = TRUE, clip = "off",
                            xlim = c(-root_age - pad_left, pad_right),
                            ylim = c(0, ape::Ntip(tree) + 1),
                            pos = "bottom", abbrv = FALSE, size = 3,
                            height = grid::unit(1, "line")) +
        ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1.5, 0.5, 0.5, "cm"))
    }
  }

  if (save) {
    x_range <- diff(range(tree_plot$data$x))

    if (layout == "circular") {
      # ggtree stores circular-layout node coordinates pre-polar-transform:
      # $data$x is still root-to-tip depth (the actual radius once rendered)
      # and $data$y is just each tip's 1..Ntip row index, not a spatial
      # length -- so only x_range (the radius) determines the rendered
      # tree's physical size; using y_range here would scale the canvas to
      # thousands of inches on trees with many tips. A circular tree radiates
      # tip labels in every direction, so the saved canvas is a square built
      # from the full diameter (2x radius) plus equal padding on all sides.
      side <- 2 * x_range + 2 * (tip_label_extent + char_width)
      base_width <- side
      base_height <- side
    } else {
      base_width <- x_range + tip_label_extent + char_width
      # Row height per tip only needs to be large enough to keep labels from
      # overlapping; with no tip labels drawn (e.g. large overview trees with
      # hundreds/thousands of tips), a much smaller per-tip factor keeps the
      # canvas a manageable size instead of scaling to hundreds of inches.
      tip_row_height <- if (show.tip.label) 15 / 100 else 3 / 100
      base_height <- max(ape::Ntip(tree) * tip_row_height, 5)
    }

    .save_tree(fulltree = tree_plot,
              fullname = fullname,
              dpi = dpi,
              format = format,
              base_height = base_height,
              base_width = base_width)
  }

  return(tree_plot)
}
