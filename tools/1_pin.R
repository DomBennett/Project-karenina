# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# LIBS
library (MoreTreeTools)
library (paleobioDB)
source (file.path ('tools', 'palaeo_tools.R'))

# INPUT
# use hominoids for test
data('hominoids')
tree <- hominoids
rm (hominoids)
tree <- multi2di (tree)

records <-  pbdb_occurrences (limit=10,
                              base_name="hominoidea", vocab="pbdb",
                              show=c("phylo", "ident"))
# remove non-species
records <- records[records$taxon_rank == 'species',]
# merge occurrence records into species records
records$taxon_no == records$taxon_no

hist (records$early_age)

res <- data.frame (name=NA, max.age=NA, min.age=NA)
lineages <- list ()
records$binomial <- paste0 (records$genus_name, '_', records$species_name)
binomials <- unique (records$binomial)
max.age <- min.age <- rep (NA, length (binomials))
for (i in 1:length (binomials)) {
  binomial <- binomials[i]
  pull <- records$binomial == binomial
  lineages[[i]] <- records[i, c ('phylum', 'class', 'order', 'family', 'genus_name', 'species_name')]
  max.age[i] <- max (records$early_age[pull])
  min.age[i] <- min (records$late_age[pull])
}
res <- data.frame (name=binomials, min.age, max.age)
.mnGetLSR (qry=lineages[[1]], sbjcts=resolve.list$lineage)

records <- data.frame (records, stringsAsFactors=FALSE)

as.character (records[1, c ('phylum', 'class', 'order', 'family', 'genus_name', 'species_name')])


pinNames <- function (tree, names, fuzzy=TRUE, datasource=4,
                      iterations=1, resolve.list=NULL) {
  # SAFETY CHECK
  if (iterations < 1) {
    stop ('Iterations must be >=1')
  }
  if (!is.vector (names) && length (names) <= 1) {
    stop ('Names must be a vector of more than 1 character string')
  }
  if (!class (tree) %in%  c ('phylo', 'multiPhylo')) {
    stop ('Tree must be phylo or multiPhylo')
  }
  #TODO update addTip to handle node.label, for now drop node.label
  # INIT
  # drop underscores, check for branch lengths and drop node labels
  tree <- .mnClean (tree)
  names <- gsub ('_', ' ', names)
  # find matching and non-matching names
  tip.labels <- .mnGetNames (tree)
  matching.names <- names[names %in% tip.labels]
  nonmatching.names <- names[!names %in% tip.labels]
  # stop now if all names match names in tree or fuzzy is false
  if (!fuzzy || length (matching.names) == length (names)) {
    return (.mnEarlyReturn (tree, names, iterations))
  }
  # VARIABLES AND ENVIRONMENTS
  # place all parameters in an environment to save arguments
  paraenv <- new.env (parent=emptyenv ())
  if (class (tree) == 'multiPhylo') {
    paraenv$trees <- tree
    paraenv$start.tree <- tree[[sample (1:length (tree), 1)]]
  } else {
    paraenv$start.tree <- tree
  }
  paraenv$grow.tree <- paraenv$start.tree
  paraenv$datasource <- datasource
  paraenv$matching.names <- matching.names
  paraenv$names <- names
  # get query and subject resolved names
  if (!is.null (resolve.list)) {
    # if resolve.list provided, no need to do any searches
    # unpack resolved names into qrylist and sbjctenv
    pull <- resolve.list$resolved$search.name %in% nonmatching.names
    if (sum (pull) < 1) {
      return (.mnEarlyReturn (tree, names, iterations))
    }
    qrylist <- list ()
    qrylist$resolved <- resolve.list$resolved[pull, ]
    qrylist$lineages <- resolve.list$lineages[pull]
    sbjctenv <- new.env (parent=emptyenv ())
    sbjctenv$resolved <- resolve.list$resolved
    sbjctenv$lineages <- resolve.list$lineages
    paraenv$deja.vues <- tree$tip.label
  } else {
    # hold query name resolution results in a list
    qrylist <- .mnResolve (names=nonmatching.names, paraenv=paraenv)
    if (!nrow (qrylist$resolved) > 0) {  # stop now if none resolved
      return (.mnEarlyReturn (tree, names, iterations))
    }
    # hold subject name resolution results in a single env
    sbjctenv <- new.env (parent=emptyenv ())
    .mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
    # add qrylist results to sbjctenv
    pull <- !as.vector (qrylist$resolved$search.name) %in%
      as.vector (sbjctenv$resolved$search.name)
    sbjctenv$resolved <- rbind (sbjctenv$resolved, qrylist$resolved[pull, ])
    sbjctenv$lineages <- c (sbjctenv$lineages, qrylist$lineages[pull])
  }
  # hold all resulting trees and stats in a single env
  resenv <- new.env (parent=emptyenv ())
  resenv$trees <- list ()
  # ITERATE
  # loop through for each iteration, write results to resenv
  m_ply (.data=data.frame(iteration=1:iterations), .fun=.mnMap,
         resenv=resenv, qrylist=qrylist, sbjctenv=sbjctenv,
         paraenv=paraenv)
  # RETURN
  if (length (resenv$trees) > 1) {
    trees <- resenv$trees
    class (trees) <- 'multiPhylo'
  } else {
    trees <- resenv$trees[[1]]
  }
  return (trees)
}






source("functions/TaxaResolve.R")
source("functions/EcoDataTools.R")
source("functions/LFITools.R")

cleanNames <- function (names) {
  # Remove name accronyms for name matching
  clean <- function (name) {
    name <- gsub ("\\s+cf\\.", "", name)
    name <- gsub ("\\s+sp\\.", "", name)
    name <- gsub ("\\s+ex\\.", "", name)
    name <- gsub ("\\s+gr\\.", "", name)
    name <- gsub ("\\s+aff\\.", "", name)
    name <- gsub ("\\s+n\\.", "", name)
    name <- gsub ("\\s+gen\\.", "", name)
    name <- gsub ("\\s+\\?", "", name)
    #name <- sub ("\\s+indet\\.", "", name)
  }
  mdply (.data = data.frame (name = names), .fun = clean)[ ,2]
}

## Input
output.dir <- "0_data"
input.dir <- file.path(output.dir, "raw")
phylo <- read.tree (file.path (input.dir, "bininda.txt"))
if (!is.binary.tree (phylo)) {
  phylo <- multi2di (phylo)
}
phylo$tip.label <- sub ("_", " ", phylo$tip.label)

## Process
cat ("Retrieving fossil records ... \n")

cat ("Cleaning names ... \n")
records$name <- cleanNames (records$name)
cat ('Labelling phylogeny ... \n')
phylo <- labelNodes (phylo)
cat ('Adding node ages ...\n')
phylo <- addNodeAges (phylo)
cat ('Adding fossils to phylogeny ...\n')
new.phylo <- addFossilsToPhylogeny (phylo, records)
cat ('Outputting ...\n')
options ("expressions" = 20000)
write.tree (new.phylo, file.path (output.dir, "mammalia_w_fossils.tre"))