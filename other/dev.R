library(treeman)

# quick dev functions
viz <- function(tree) {
  plot(as(tree, 'phylo'))
  ape::axisPhylo()
}
nTr <- function(n) {
  tree <- ape::compute.brlen(ape::rtree(n))
  as(tree, 'TreeMan')
}

# dev functions

tree <- nTr(50)
viz(tree)