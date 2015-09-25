# D.J. Bennett
# 17/07/2015
# Tools for ED and fossils

# LIBS
source (file.path ('tools', 'pin_names.R'))

# PIN TOOLS -- requires parameters.Rd
pin <- function (tree=tree, names=binomials, lineages=lineages,
                 min.ages=min.age, max.ages=max.age,
                 iterations=iterations, resolve.list=resolve.list) {
  if (ncpus < 2) {
    pin.res <- pinNames (tree=tree, names=binomials, lineages=lineages,
                         min.ages=min.age, max.ages=max.age,
                         iterations=iterations, resolve.list=resolve.list)
    return (pin.res)
  }
  # UNIX-only parrelisation
  require (foreach)
  require (doMC)
  registerDoMC (ncpus)
  pin.res <- foreach (i=1:iterations) %dopar% {
    cat ('\n........ [', i, ']', sep='')
    pin.res <- pinNames (tree=tree, names=binomials, lineages=lineages,
                         min.ages=min.age, max.ages=max.age,
                         iterations=1, resolve.list=resolve.list)
    pin.res
  }
  class (pin.res) <- 'multiPhylo'
  pin.res
}

# SLICE TOOLS

addTips2Nodes <- function (tree) {
  # Return a tree with new tips at every node
  # All tips labelled N1-NNode
  nodedict <- getNodedict (tree)
  for (i in 1:length (nodedict)) {
    # ignore tips with only one descendant
    if (length (nodedict[[i]]) < 2) {
      next
    }
    children <- nodedict[[i]]
    node <- getParent (tree, tip=children)
    edge <- which (tree$edge[ ,1] == node)[1]
    node.age <- getAge (tree, node=node)[,2]
    tree <- addTip (tree=tree, edge=edge,
                    node.age=node.age, tip.age=node.age,
                    tip.name=names (nodedict)[i])
  }
  tree$tip.label <- paste0 ('n', 1:getSize (tree))
  tree
}

getNodedict <- function (tree) {
  # Return a list of each node and its children
  nodedict <- list ()
  nnode <- getSize (tree) + tree$Nnode
  for (i in 1:nnode) {
    node.name <- paste0 ('n', i)
    nodedict[node.name] <- list (getChildren (tree, i))
  }
  nodedict
}

calcEDBySlice <- function (tree, time.cuts) {
  # Return ED values for clades at different time slices
  # for time callibrated tree
  age <- max (diag (vcv.phylo (tree)))
  time.cuts <- time.cuts[time.cuts < age]
  all.node.labels <- paste0 ('n', 1:(length (tree$tip.label) + tree$Nnode))
  res <- matrix (nrow=length (time.cuts), ncol=length (all.node.labels))
  rownames (res) <- as.character (time.cuts)
  colnames (res) <- all.node.labels
  for (i in 1:length (time.cuts)) {
    # slice tree at interval
    sliced <- getTimeslice (tree=tree, time.slice=time.cuts[i],
                            all.node.labels=all.node.labels)
    # get ed vals
    ed.res <- calcED (sliced)
    # normalise by PD
    ed.res$ED <- ed.res$ED/sum (sliced$edge.length)
    # add to res
    indexes <- match (rownames(ed.res), all.node.labels)
    res[i,indexes] <- ed.res[,1]
  }
  return (res)
}


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