# Sourced parameters
ncpus <- time_cuts <- overwrite <- iterations <- parent <- treefile <- NULL

# PIN TOOLS
assembleRecords <- function(records, taxonomy, tree_tips, ...) {
  # parse names, linages and ages for pinning
  lineages <- list()
  records$binomial <- paste0(records$genus_name, ' ', records$species_name)
  records$binomial <- gsub('\\s', '_', records$binomial)
  records$binomial <- gsub('[^a-zA-Z0-9_]', '', records$binomial)
  records <- records[!records$binomial %in% tree_tips, ]
  binomials <- unique(records$binomial)
  max.age <- min.age <- rep(NA, length(binomials))
  .calc <- function(i) {
    binomial <- binomials[i]
    pull <- records$binomial == binomial
    # get age range for species
    max.age[i] <- max(records$early_age[pull])
    min.age[i] <- min(records$late_age[pull])
    # extract PDBD taxonomy
    lineage <- records[which(pull)[1], taxonomy, drop = TRUE]
    lineage <- gsub('\\s', '_', as.vector(unlist(lineage)))
    lineage <- gsub('[^a-zA-Z0-9_]', '', lineage)
    lineage <- lineage[!duplicated(lineage)]
    lineage <- lineage[!is.na(lineage)]
    lineages[[i]] <- lineage[-(length(lineage))]  # rm species name
    # push to parent environment
    max.age <<- max.age
    min.age <<- min.age
    lineages <<- lineages
  }
  l.data <- data.frame(i = 1:length(binomials))
  plyr::m_ply(.data = l.data, .fun = .calc, ...)
  list('lineages' = lineages, 'max.age' = max.age, 'min.age' = min.age,
        'binomials' = binomials)
}

pinParallel <- function(tree, tids, lngs, min_ages, max_ages, pinfolder,
                        tree_age) {
  if (!file.exists(pinfolder)) {
    dir.create(pinfolder)
  }
  if (ncpus < 2) {
    end_ages <- sapply(1:length(tids), function(x){
      runif(1, min = min_ages[x], max = max_ages[x])
    })
    tree <- pinTips(tree, tids, lngs, end_ages, tree_age)
    save(tree, file = file.path(pinfolder, '1.RData'))
    return(NULL)
  }
  # UNIX-only parrelisation
  require(foreach)
  require(doMC)
  registerDoMC(ncpus)
  outfiles <- file.path(pinfolder, paste0(iterations, '.RData'))
  if (!overwrite) {
    iterations <- iterations[!file.exists(outfiles)]
    outfiles <- outfiles[!file.exists(outfiles)]
  }
  if (length(outfiles) == 0) {
    message('All iterations already run.')
    return(invisible(NULL))
  }
  i <- NULL
  res <- foreach(i = 1:length(iterations)) %dopar% {
    cat('\n........ [', iterations[i], ']', sep = '')
    end_ages <- sapply(1:length(tids), function(x){
      runif(1, min = min_ages[x], max = max_ages[x])
    })
    # not necessarily efficient to run in parallel
    # http://stackoverflow.com/questions/22104858/is-it-a-good-idea-to-read-write-files-in-parallel
    # .... but better to save progress
    # re-jig order in case it impacts pinning
    tree <- pinTips(tree, tids, lngs, end_ages, tree_age)
    save(tree, file = outfiles[i])
  }
  invisible(res)
}
