
library(tidygraph)
library(ggraph)
library(ggplot2)

# Create edge list
edges <- data.frame(
  from = c("A", "B", "C", "D", "E", "F", "G", "G", "I", "J", "K"),
  to = c("C", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
)

# Create node list with labels
nodes <- data.frame(
  name = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
  label = c(
    "mineSeq()\nSequence retrieval",
    "minePlastome() /\nmineMitochondrion()\nTargeted locus mining",
    "combineFASTA()\nCombine FASTA files",
    "alignSeqs()\nMultiple sequence alignment",
    "convertAlign()\nAlignment conversion",
    "ape::read.nexus.data()\nAlignment import",
    "catfullGenes() /\ncatmultGenes()\nConcatenation",
    "dropSeq()\nRemove redundant\naccessions",
    "writeNexus() /\nwritePhylip()\nExport datasets",
    "evomodelTest()\nModel selection and\nMrBayes block",
    "mrbayesRun()\nBayesian inference",
    "plotPhylo()\nTree visualization"
  ),
  stage = c(rep("retrieval", 2), rep("preparation", 4),
            rep("concatenation", 3), rep("analysis", 3))
)

# Create graph
graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

# Plot with ggraph
p <- ggraph(graph, layout = "sugiyama") +
  geom_edge_diagonal(aes(color = "workflow"),
                     arrow = arrow(length = unit(4, 'mm')),
                     end_cap = circle(3, 'mm'),
                     start_cap = circle(3, 'mm')) +
  geom_node_point(aes(color = stage), size = 20) +
  geom_node_text(aes(label = label), size = 3) +
  scale_color_manual(values = c(
    "retrieval" = "#E8F4FD",
    "preparation" = "#FFF2CC",
    "concatenation" = "#D5E8D4",
    "analysis" = "#E1D5E7"
  )) +
  theme_graph() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "catGenes Workflow Overview")

# Display
print(p)

# Save
ggsave("inst/figures/catGenes_workflow_ggraph.png", p, width = 14, height = 10, dpi = 300)
