# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# PARAMETERS
parent <- "hominoidea"  # name of parent clade for species in tree
iterations <- 100  # number of trees

# DIRS
output.dir <- '1_pin'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# LIBS
library (MoreTreeTools)
library (paleobioDB)
source (file.path ('tools', 'pin_names.R'))

# INPUT
# use hominoids for test
data('hominoids')
tree <- hominoids
rm (hominoids)
tree <- multi2di (tree)
tree <- drop.tip (tree, tip='Macaca mulatta')

# PALEODB
records <-  pbdb_occurrences (limit='all',
                              base_name=parent, vocab="pbdb",
                              show=c("phylo", "ident"))

# EXTRACT NAMES + LINEAGES
taxonomy <- c ('phylum', 'class', 'order', 'family', 'genus_name', 'species_name')
lineages <- list ()
records$binomial <- paste0 (records$genus_name, '_', records$species_name)
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
write.csv (records, file=file.path (output.dir, 'palaeo_records.csv'),
           row.names=FALSE)
write.tree (pin.res, file=file.path (output.dir, 'pinned.tre'))
