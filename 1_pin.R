# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# TIMESTAMP
cat (paste0 ('\nPinning started at [', Sys.time (), ']'))

# PARAMETERS
load ('parameters.Rd')

# DIRS
input.dir <- '0_data'
output.dir <- '1_pin'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# LIBS
library (MoreTreeTools)
library (paleobioDB)
source (file.path ('tools', 'palaeo_tools.R'))

# INPUT
tree <- read.tree (file.path (input.dir, treefile))

# PALEODB
cat ('\nRetrieving records ....')
fossilfile <- paste0 (parent, '_records.csv')
if (overwrite | !file.exists (file.path (output.dir, fossilfile))) {
  # get
  cat ('\n.... searching PBDB')
  records <-  pbdb_occurrences (limit='all',
                                base_name=parent, vocab="pbdb",
                                show=c("phylo", "ident"))
  # write out
  write.csv (records, file=file.path (output.dir, fossilfile),
             row.names=FALSE)
} else {
  cat ('\n.... pre-loading file')
  records <- read.csv (file=file.path (output.dir, fossilfile))
}
cat ('\nDone.')

# EXTRACT NAMES + LINEAGES
cat ('\nAssembling records for pinning ....')
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
cat ('\nDone.')

# PIN
cat ('\nPinning ....')
names <- tree$tip.label
rlfile <- paste0 (parent, '_resolvedlist.Rd')
if (overwrite | !file.exists (file.path (output.dir, rlfile))) {
  resolve.list <- mapNamesPreDownload(names, datasource = 4)
  save (resolve.list, file=file.path (output.dir, rlfile))
} else {
  load (file=file.path (output.dir, rlfile))
}
cat ('\n.... real')
pin.res <- pin (tree=tree, names=binomials, lineages=lineages,
                min.ages=min.age, max.ages=max.age,
                iterations=iterations, resolve.list=resolve.list)
pinfile <- paste0 (parent, '_pinned.tre')
write.tree (pin.res, file=file.path (output.dir, pinfile))
# randomise ages
randis <- sample (1:length (lineages))
min.age <- min.age[randis]
max.age <- max.age[randis]
# randomise linages
lineages <- lineages[sample (1:length (lineages))]
cat ('\n.... random')
# rand.res <- pin (tree=tree, names=binomials, lineages=lineages,
#                  min.ages=min.age, max.ages=max.age,
#                  iterations=iterations, resolve.list=resolve.list)
# randfile <- paste0 (parent, '_rand.tre')
# write.tree (rand.res, file=file.path (output.dir, randfile))
cat ('\nDone.')

# TIMESTAMP
cat (paste0 ('\nPinning finished at [', Sys.time (), ']'))