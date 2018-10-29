# Sourced parameters
ncpus <- time_cuts <- overwrite <- iterations <- parent <- treefile <- NULL
source(file.path('tools', 'pin_tips.R'))

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
  # UNIX-only parellelisation
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
  ntips <- foreach(i = 1:length(iterations)) %dopar% {
    cat('\n........ [', iterations[i], ']', sep = '')
    # end_ages <- sapply(1:length(tids), function(x){
    #   runif(1, min = min_ages[x], max = max_ages[x])
    # })
    # not necessarily efficient to run in parallel
    # http://stackoverflow.com/questions/22104858/is-it-a-good-idea-to-read-write-files-in-parallel
    # .... but better to save progress
    # re-jig order in case it impacts pinning
    tree <- suppressMessages(pinTips2(tree, tids, lngs, min_ages, tree_age))
    save(tree, file = outfiles[i])
    tree@ntips
  }
  ntips
}

pinTips_rand <- function(tree, tids, end_ages, tree_age) {
  spn_data <- matrix(NA, nrow = (length(tids) * 2) + tree@nall, ncol = 2)
  colnames(spn_data) <- c("start", "end")
  tmp_spn_data <- getSpnsAge(tree, tree@all, tree_age)
  rownames(spn_data) <- c(tree@all, tids, paste0("p_", tids))
  spn_data[tree@all, "start"] <- tmp_spn_data[["start"]]
  spn_data[tree@all, "end"] <- tmp_spn_data[["end"]]
  rm(tmp_spn_data)
  for (i in seq_along(tids)) {
    tid <- tids[[i]]
    end_age <- end_ages[[i]]
    pull <- unname(which(spn_data[ ,'start'] > end_age))
    probs <- unname(spn_data[pull, 'start'] - spn_data[pull, 'end'])
    randi <- sample(pull, 1, prob = probs/max(probs))
    sid <- rownames(spn_data)[randi]
    pssbl_ages <- getSpnAge(tree, sid, tree_age)
    if (end_age > pssbl_ages[['end']]) {
      strt_age <- runif(1, end_age, pssbl_ages[['start']])
    } else {
      strt_age <- runif(1, pssbl_ages[['end']], pssbl_ages[['start']])
    }
    success <- tryCatch(expr = {
      tree <- addTip(tree = tree, tid = tid, sid = sid, strt_age = strt_age,
                     end_age = end_age, tree_age = tree_age)
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (!success) {
      next
    }
    # update spn_data with new info
    pid <- paste0("p_", tid, sep = "")
    tid_spn <- getSpnAge(tree, tid, tree_age)
    spn_data[tid, "start"] <- tid_spn[, "start"]
    spn_data[tid, "end"] <- tid_spn[, "end"]
    pid_spn <- getSpnAge(tree, pid, tree_age)
    spn_data[pid, "start"] <- pid_spn[, "start"]
    spn_data[pid, "end"] <- pid_spn[, "end"]
    sid_spn <- getSpnAge(tree, sid, tree_age)
    spn_data[sid, "start"] <- sid_spn[, "start"]
    spn_data[sid, "end"] <- sid_spn[, "end"]
  }
  tree
}

pinParallelRand <- function(tree, tids, pinfolder, min_ages, tree_age) {
  if (!file.exists(pinfolder)) {
    dir.create(pinfolder)
  }
  outfiles <- file.path(pinfolder, paste0(iterations, '.RData'))
  if (!overwrite) {
    iterations <- iterations[!file.exists(outfiles)]
    outfiles <- outfiles[!file.exists(outfiles)]
  }
  if (length(outfiles) == 0) {
    message('All iterations already run.')
    return(invisible(NULL))
  }
  ntips <- foreach(i = 1:length(iterations)) %dopar% {
    cat('\n........ [', iterations[i], ']', sep = '')
    tree <- pinTips_rand(tree, tids, min_ages, tree_age)
    save(tree, file = outfiles[i])
    tree@ntips
  }
  ntips
}
