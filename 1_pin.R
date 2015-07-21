# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# PARAMETERS
parent <- "mammalia"  # name of parent clade for species in tree
iterations <- 100  # number of trees
overwrite <- FALSE

# DIRS
input.dir <- '0_data'
output.dir <- '1_pin'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# LIBS
library (MoreTreeTools)
library (paleobioDB)
source (file.path ('tools', 'pin_names.R'))

# INPUT
tree <- read.tree (file.path (input.dir, 'bininda_mammalia.tre'))
tree <- multi2di (tree)

# PALEODB
fossilfile <- paste0 (parent, '_records.csv')
if (overwrite | !file.exists (file.path (output.dir, fossilfile))) {
  records <-  pbdb_occurrences (limit='all',
                                base_name=parent, vocab="pbdb",
                                show=c("phylo", "ident"))
} else {
  records <- read.csv (file=file.path (output.dir, fossilfile))
}

# EXTRACT NAMES + LINEAGES
taxonomy <- c ('phylum', 'class', 'order', 'family', 'genus_name', 'species_name')
lineages <- list ()
records$binomial <- paste0 (records$genus_name, ' ', records$species_name)
records <- records[!records$binomial %in% tree$tip.label, ]
binomials <- unique (records$binomial)
max.age <- min.age <- rep (NA, length (binomials))
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

# PIN
pin.res <- pinNames (tree=tree, names=binomials, lineages=lineages,
                     min.ages=min.age, max.ages=max.age,
                     iterations=iterations)

# OUTPUT
write.csv (records, file=file.path (output.dir, fossilfile),
           row.names=FALSE)
write.tree (pin.res, file=file.path (output.dir, 'pinned.tre'))
