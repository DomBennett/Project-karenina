# PIN TOOLS
assembleRecords <- function (records, taxonomy, ...) {
  # parse names, linages and ages for pinning
  lineages <- list ()
  records$binomial <- paste0 (records$genus_name, ' ', records$species_name)
  records <- records[!records$binomial %in% tree_tips, ]
  binomials <- unique (records$binomial)
  max.age <- min.age <- rep (NA, length (binomials))
  .calc <- function (i) {
    binomial <- binomials[i]
    pull <- records$binomial == binomial
    # get age range for species
    max.age[i] <- max(records$early_age[pull])
    min.age[i] <- min(records$late_age[pull])
    # extract PDBD taxonomy
    lineage <- records[which (pull)[1], taxonomy, drop=TRUE]
    lineages[[i]] <- as.vector (unlist (lineage))
    # push to parent environment
    max.age <<- max.age
    min.age <<- min.age
    lineages <<- lineages
  }
  l.data <- data.frame (i=1:length (binomials))
  plyr::m_ply (.data=l.data, .fun=.calc, ...)
  list ('lineages'=lineages, 'max.age'=max.age, 'min.age'=min.age,
        'binomials'=binomials)
}

pinParallel <- function (tree, tids, lngs, min_ages, max_ages, outfile) {
  if(file.exists(outfile)) {
    file.remove(outfile)
  }
  if (ncpus < 2) {
    end_ages <- sapply(1:length(tids), function(x){
      runif(1, min=min_ages[x], max=max_ages[x])
    })
    tree <- pinTips(tree, tids, lngs, end_ages, tree_age)
    writeTree(tree, file=outfile)
  }
  # UNIX-only parrelisation
  require (foreach)
  require (doMC)
  registerDoMC (ncpus)
  trees <- foreach (i=1:iterations) %dopar% {
    cat ('\n........ [', i, ']', sep='')
    end_ages <- sapply(1:length(tids), function(x){
      runif(1, min=min_ages[x], max=max_ages[x])
    })
    pinTips(tree, tids, lngs, end_ages, tree_age)
  }
  for(tree in trees) {
    # not efficient to run in parallel
    # http://stackoverflow.com/questions/22104858/is-it-a-good-idea-to-read-write-files-in-parallel
    writeTree(tree, file=outfile, append=TRUE)
  }
}
