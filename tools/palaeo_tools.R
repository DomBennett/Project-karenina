## No copyright, no warranty
## Dominic John Bennett
## Functions for calculating LFI
## 07/02/2014

## Dependencies
require (ape)
require (geiger)
require (plyr)
library (graphics)

## Functions
getTimeslice <- function (tree, time.slice,
                          all.node.ages=getAge (tree)$age,
                          all.node.labels=paste0 ('c', 1:(length (tree$tip.label) +
                                                            tree$Nnode))) {
  # identify all nodes and edges that need to be kept
  # all edges that have 1 node on older than time.slice
  keep.nodes <- which (all.node.ages >= time.slice)
  keep.edges <- which (tree$edge[ ,1] %in% keep.nodes)
  # remove any tip nodes that don't pass through timeslice
  tip.edges <- keep.edges [tree$edge[keep.edges, 2] <= length (tree$tip.label)]
  tip.node.ages <- all.node.ages[tree$edge[tip.edges, 2]]
  drop.tip.edges <- tip.edges[tip.node.ages >= time.slice]
  keep.edges <- keep.edges[!keep.edges %in% drop.tip.edges]
  # create an edge matrix
  edge.matrix <- tree$edge[keep.edges, ]
  edge.lengths <- tree$edge.length[keep.edges]
  tips <- edge.matrix[!edge.matrix[ ,2] %in% edge.matrix[ ,1],2]
  tip.labels <- all.node.labels[tips]
  n.node <- as.integer (length (unique (edge.matrix[ ,1])))
  # relabel all nodes
  new.edge.matrix <- edge.matrix
  new.edge.matrix[edge.matrix[ ,2] %in% tips,2] <- 1:length (tips)
  new.edge.matrix[!new.edge.matrix[ ,1] %in% new.edge.matrix[ ,2],1] <- length (tips) + 1
  old.nodes <- unique (new.edge.matrix[new.edge.matrix[ ,1] > length (tips) + 1,1])
  new.nodes <- (length (tips) + 2):(length (old.nodes) + (length (tips) + 1))
  for (i in 1:length (old.nodes)) {
    new.edge.matrix[new.edge.matrix[ ,1] == old.nodes[i],1] <- new.nodes[i]
    new.edge.matrix[new.edge.matrix[ ,2] == old.nodes[i],2] <- new.nodes[i]
  }
  storage.mode (new.edge.matrix) <- "integer"
  # create new tree object
  new.tree <- list (edge = new.edge.matrix, tip.label = tip.labels,
                     edge.length = edge.lengths, Nnode = n.node)
  class (new.tree) <- 'phylo'
  # identify corresponding nodes
  old.nodes <- unique (c (edge.matrix[ ,1], edge.matrix[ ,2]))
  new.nodes <- unique (c (new.edge.matrix[ ,1], new.edge.matrix[ ,2]))
  all.node.ages <- all.node.ages[old.nodes]
  new.all.node.ages <- all.node.ages[order (new.nodes)]
  # identify all tip nodes and re-lengthen
  for (i in 1:length (new.tree$tip.label)) {
    diff <- time.slice - new.all.node.ages[i]
    edge <- which (new.tree$edge[ ,2] == i)
    edge.length <- new.tree$edge.length[edge]
    new.edge.length <- edge.length - diff
    new.tree$edge.length[edge] <- new.edge.length
  }
  collapse.singles (new.tree)
}