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
                                show=c("phylo", "ident"),
                                stringsAsFactors=FALSE)
  # write out
  write.csv (records, file=file.path (output.dir, fossilfile),
             row.names=FALSE)
} else {
  cat ('\n.... pre-loading file')
  records <- read.csv (file=file.path (output.dir, fossilfile),
                       stringsAsFactors=FALSE)
}
cat ('\nDone.')

# EXTRACT NAMES + LINEAGES
cat ('\nAssembling records for pinning ....')
# specify ranks in linages for pinning
taxonomy <- c ('phylum', 'class', 'order', 'family', 'genus_name', 'species_name')
records.obj <- assembleRecords (records, taxonomy, .progress='time')
lineages <- records.obj[['lineages']]
length (lineages)
max.age <- records.obj[['max.age']]
min.age <- records.obj[['min.age']]
binomials <- records.obj[['binomials']]
length (binomials)
rm (records.obj)

fossils <- list("binomial"=binomials, "min_age"=min.age,
                "max_age"=max.age, "lineage"=lineages)
tree_taxonomy <- resolve.list

min.age
max.age
lineages
resolve.list

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
pinfile <- file.path (output.dir, paste0 (parent, '_pinned.tre'))
pin.res <- pin (tree=tree, names=binomials, lineages=lineages,
                min.ages=min.age, max.ages=max.age,
                iterations=iterations, resolve.list=resolve.list,
                outfile=pinfile)
# randomise ages
# randis <- sample (1:length (lineages))
# min.age <- min.age[randis]
# max.age <- max.age[randis]
# # randomise linages
# lineages <- lineages[sample (1:length (lineages))]
# cat ('\n.... random')
# rand.res <- pin (tree=tree, names=binomials, lineages=lineages,
#                  min.ages=min.age, max.ages=max.age,
#                  iterations=iterations, resolve.list=resolve.list)
# randfile <- paste0 (parent, '_rand.tre')
# write.tree (rand.res, file=file.path (output.dir, randfile))
cat ('\nDone.')

# TIMESTAMP
cat (paste0 ('\nPinning finished at [', Sys.time (), ']'))