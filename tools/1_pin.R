# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# LIBS
library (MoreTreeTools)
library (paleobioDB)]
source (file.path ('tools', 'palaeo_tools.R'))

calcMeanDiffs <- function (ed.changes) {
  # calculate the mean differences in relative ED
  difference <- starting <- rep (NA, ncol (ed.changes))
  for (i in 1:ncol (ed.changes)) {
    starting[i] <- mean (ed.changes[ ,i], na.rm=TRUE)
    difference[i] <- mean (ed.changes[-nrow(ed.changes),i] - ed.changes[-1,i], na.rm=TRUE)
  }
  res <- data.frame (mean.start=starting, mean.difference=difference)
  rownames (res) <- colnames (ed.changes)
  res
}

# INPUT
# use hominoids for test
data('hominoids')
tree <- hominoids
rm (hominoids)
tree <- multi2di (tree)
tree <- drop.tip (tree, tip='Macaca mulatta')
# get all ape records
records <-  pbdb_occurrences (limit='all',
                              base_name="hominoidea", vocab="pbdb",
                              show=c("phylo", "ident"))
lineages <- list ()
records$binomial <- paste0 (records$genus_name, '_', records$species_name)
binomials <- unique (records$binomial)
max.age <- min.age <- rep (NA, length (binomials))
taxonomy <- c ('phylum', 'class', 'order', 'family', 'genus_name', 'species_name')
for (i in 1:length (binomials)) {
  binomial <- binomials[i]
  pull <- records$binomial == binomial
  # get age range for species
  max.age[i] <- max (records$early_age[pull])
  min.age[i] <- min (records$late_age[pull])
  # extract PDBD taxonomy
  lineage <- records[which(pull)[1], taxonomy, drop=TRUE]
  lineages[[i]] <- as.vector(unlist (lineage))
}

pin.res <- pinNames (tree=tree, names=binomials, lineages=lineages,
                     min.ages=min.age, max.ages=max.age)
plot (pin.res)
# take timeslices of pin.res
steps <- 10
age <- max (diag (vcv.phylo (pin.res)))
interval <- age/steps
intervals <- seq (from=interval, to=age, by=interval)
all.node.labels <- paste0 ('n', 1:(length (pin.res$tip.label) + pin.res$Nnode))
res <- matrix (nrow=steps, ncol=length (all.node.labels))
colnames (res) <- all.node.labels
for (i in 1:steps) {
  # slice tree at interval
  sliced <- getTimeslice (tree=pin.res, time.slice=intervals[i],
                          all.node.labels=all.node.labels)
  # get ed vals
  ed.res <- calcED (sliced)
  # normalise by PD
  ed.res$ED <- ed.res$ED/sum (sliced$edge.length)
  # add to res
  indexes <- match (rownames(ed.res), all.node.labels)
  res[i,indexes] <- ed.res[,1]
}
# how does the starting ED predict the future ED?
diff.res <- calcMeanDiffs (res)
plot (abs(diff.res$mean.difference) ~ diff.res$mean.start)


# Leftover

plot (pin.res)

pin.res$all.node.label <- 1:(length (pin.res$tip.label) + pin.res$Nnode)
par(mfrow=c(1,2))
display.tree <- pin.res
display.tree$tip.label <- display.tree$all.node.label[1:length(display.tree$tip.label)]
plot(display.tree);nodelabels(text=display.tree$all.node.label[
  (length(display.tree$tip.label)+1):length(display.tree$all.node.label)])
axisPhylo()
# does starting ED predict change in ED?
ed.diffs <- ed.res$t1 - ed.res$t2
plot (ed.diffs ~ ed.res$t2)

clade.refs <- list ()
all.node.labels <- paste0 ('c', 1:(length (tree$tip.label) + tree$Nnode))
clade.refs['c1'] <- list('', children=getChildren(pin.res))

plot (sliced);axisPhylo()




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