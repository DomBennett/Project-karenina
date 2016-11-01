library(treeman)

# Load tree
load("~/Desktop/mammal.RData")
tips(tree) <- gsub("_", " ", tips(tree))

# 1. taxonomise
#  - search taxonomy for all names in tree
#  - get children by node
#  - find highest shared taxonomic group per node based on taxonomies of children
#  - assign to taxonym

.findClade <- function(lineages) {
  # for a list of lineages, find the lowest clade shared by all
  .mtch <- function(qry) {
    sbj <<- sbj[sbj %in% qry]
    NULL
  }
  lineages <- lineages[!is.na(lineages)]
  sbj <- lineages[[1]]
  lapply(lineages, .mtch)
  sbj[length(sbj)]
}

.addTxnym <- function(nd, lineages) {
  if(length(nd$children) > 0) {
    mtchs <- match(nd$children, resolved$search.name)
    txnym <- .findClade(lineages[mtchs])
    nd$taxonym <- txnym
  } else {
    nd$taxonym <- strsplit(nd$id, "\\s")[[1]][1]
  }
  nd
}
.rmvWldCrds <- function(l) {
  if(!parent %in% l) {
    l <- NA
  }
  l
}
parent <- "Mammalia"
lineages <- lapply(lineages, .rmvWldCrds)
tree@nodelist <- lapply(tree@nodelist, .addTxnym, lineage=lineages)
tree@nodelist[[sample(1:length(tree@nodelist), 1)]]  # check random nodes
# TODO: fix error with the postnodes, there are too many of them sometimes + NAs. Must be to do with renaming.

load(file="1_pin/mammalia_resolvedlist.Rd")
resolved <- resolve.list$resolved
lineages <- resolve.list$lineages





# 2. download fossil records
#  - as is done
# 3. pin fossils on to tree

children
ages <- getNodesAge(tree)
ages[1:10,]
tree[['Xenomys_nelsoni']]
