# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# TIMESTAMP
cat (paste0 ('\nPinning started at [', Sys.time (), ']'))

# PARAMETERS
source('parameters.R')

# DIRS
input.dir <- '0_data'
output.dir <- '1_pin'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# LIBS
library(treeman)
library(paleobioDB)
source(file.path ('tools', 'pin_tools.R'))

# READ TREE
cat ('\nResearching tree ....')
tree <- readTree(file.path (input.dir, treefile),
                 wndmtrx=FALSE)
tree_age <- getAge(tree)
tree_tips <- tree['tips']
txnyms <- searchTxnyms(tree, cache=TRUE)
tree <- setTxnyms(tree, txnyms)
cat ('\nDone.')

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
records.obj <- assembleRecords(records, taxonomy, .progress='time')
lineages <- records.obj[['lineages']]
max_age <- records.obj[['max.age']]
min_age <- records.obj[['min.age']]
binomials <- records.obj[['binomials']]
rm (records.obj)
cat('Done. Discovered [', length (binomials), '] records.')

# PIN
cat ('\nPinning ....')
cat ('\n.... real data')
pinfile <- file.path (output.dir, paste0 (parent, '_real.tre'))
pinParallel(tree, tids=binomials, lngs=lineages, min_ages=min_age,
            max_ages=max_age, outfile=pinfile)
cat ('\n.... random lineages')
pinfile <- file.path (output.dir, paste0 (parent, '_rndm_lngs.tre'))
pinParallel(tree, tids=binomials, lngs=sample(lineages),
            min_ages=min_age, max_ages=max_age, outfile=pinfile)
cat ('\n.... random ages')
randis <- sample(1:length (lineages))
pinfile <- file.path (output.dir, paste0 (parent, '_rndm_ages.tre'))
pinParallel(tree, tids=binomials, lngs=lineages,
            min_ages=min_age[randis],
            max_ages=max_age[randis],
            outfile=pinfile)
cat ('\nDone.')

# TIMESTAMP
cat (paste0 ('\nPinning finished at [', Sys.time (), ']'))