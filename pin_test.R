# Testing pinning using mammals

# TIMESTAMP
cat (paste0 ('\nPinning started at [', Sys.time (), ']'))

# DIRS
input.dir <- '0_data'
output.dir <- '1_pin'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# LIBS
library(treeman)

# FUNCTIONS
assembleRecords <- function (records, taxonomy, ...) {
  # parse names, linages and ages for pinning
  lineages <- list ()
  records$binomial <- paste0 (records$genus_name, '_', records$species_name)
  records <- records[!records$binomial %in% tree['tips'], ]
  binomials <- unique (records$binomial)
  max.age <- min.age <- rep (NA, length (binomials))
  .calc <- function (i) {
    binomial <- binomials[i]
    pull <- records$binomial == binomial
    # get age range for species
    max.age[i] <- max (records$early_age[pull])
    min.age[i] <- min (records$late_age[pull])
    # extract PDBD taxonomy
    lineage <- records[which (pull)[1], taxonomy, drop=TRUE]
    lineages[[i]] <- as.vector (unlist (lineage))
    # push to parent environment
    max.age <<- max.age
    min.age <<- min.age
    lineages <<- lineages
  }
  l.data <- data.frame (i=1:length (binomials))
  m_ply (.data=l.data, .fun=.calc, ...)
  list ('lineages'=lineages, 'max.age'=max.age, 'min.age'=min.age,
        'binomials'=binomials)
}


pin <- function (tree=tree, tids, lngs,
                 min_ages, max_ages, iterations,
                 outfolder) {
  .run <- function(i) {
    cat ('\n........ [', i, ']', sep='')
    # randomise order to prevent order bias
    rndi <- sample(1:length(tids))
    outfile <- file.path(outfolder, paste0('i', i, '.RData'))
    ends <- sapply(1:length(min_ages),
                   function(i) runif(min=min_ages[i], max=max_ages[i], n=1))
    pinned <- pinTips(tree, tids=tids[rndi], lngs=lngs[rndi], ends=ends[rndi])
    save(pinned, file=outfile)
  }
  if (ncpus < 2) {
    for(i in 1:iterations) {
      .run(i)
    }
  } else {
    # UNIX-only parrelisation
    require (foreach)
    require (doMC)
    registerDoMC (ncpus)
    foreach (i=1:iterations) %dopar% {
      .run(i)
    }
  }
}

# INPUT
data(mammals)
tree <- mammals
rm(mammals)

# PARAMETERS
ncpus <- 1
parent <- 'mammalia'
iterations <- 2

# PALEODB
cat ('\nRetrieving records ....')
fossilfile <- paste0(parent, '_records.csv')
cat ('\n.... pre-loading file')
records <- read.csv(file=file.path (output.dir, fossilfile),
                    stringsAsFactors=FALSE)
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

# PIN
outfolder <- file.path(output.dir, parent)
if (!file.exists(outfolder)) {
  dir.create(outfolder)
}
pin(tree, tids=binomials, lngs=lineages, min_ages=min.age,
    max_ages=max.age, iterations=iterations, outfolder=outfolder)
