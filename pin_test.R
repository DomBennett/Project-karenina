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
assembleRecords <- function (records, taxonomy, tips, ...) {
  # parse names, linages and ages for pinning
  lineages <- list ()
  records$binomial <- paste0 (records$genus_name, '_', records$species_name)
  records <- records[!records$binomial %in% tips, ]
  binomials <- unique (records$binomial)
  max.age <- min.age <- rep (NA, length (binomials))
  .calc <- function (i) {
    binomial <- binomials[i]
    pull <- records$binomial == binomial
    # get age range for species
    max.age[i] <- max (records$early_age[pull])
    min.age[i] <- min (records$late_age[pull])
    # extract PDBD taxonomy
    lineage <- records[which(pull)[1], taxonomy, drop=TRUE]
    lineage <- as.vector(unlist(lineage))
    lineages[[i]] <- lineage[!is.na(lineage)]
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


pin <- function (tree=tree, tids, lngs,
                 min_ages, max_ages, iterations,
                 outfolder) {
  .run <- function(i) {
    cat ('\n........ [', i, ']', sep='')
    # randomise order to prevent order bias
    outfile <- file.path(outfolder, paste0('i', i, '.RData'))
    rndi <- sample(1:length(tids))
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

pinTips <- function(tree, tids, lngs, ends) {
  .pin <- function(i) {
    # unpack
    tid <- tids[i]
    print(tid)
    lng <- lngs[[i]]
    end <- ends[i]
    for(j in length(lng):1) {
      spns <- names(txnyms)[which(txnyms %in% lng[j])]
      spns <- spns[spns != tree@root]
      if(length(spns) == 0) {
        next
      }
      spns <- c(spns, unlist(sapply(spns, function(n) tree@ndlist[[n]][['ptid']])))
      rngs <- getSpnsAge(tree, ids=spns)
      bool <- rngs[ ,'start'] > end
      if(any(bool)) {
        rngs <- rngs[bool, ]
        rngs[rngs[ ,'end'] <= end, "end"] <- end
        # pinning is based on branch length
        prbs <- rngs$start - rngs$end
        e <- as.vector(sample(rngs$spn, prob=prbs, size=1))
        e_i <- which(rngs$spn == e)
        start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
        if(j != length(lng)) {
          tip_txnym <- lng[j+1]
        } else {
          tip_txnym <- lng[j]
        }
        pid <- paste0('p_', tid, sep='')
        tree <- addTip(tree, tid=tid, sid=e, start=start, end=end,
                       pid=pid)
        tree@ndlst[[tid]][['txnym']] <- tip_txnym
        tree@ndlst[[pid]][['txnym']] <- lng[j]
        # add to txnyms list
        txnyms[[tid]] <<- tip_txnym
        # push out
        tree <<- tree
        break
      }
    }
  }
  .getTxnyms <- function(txnym, ...) {
    txnym
  }
  txnyms <- plyr::mlply(tree@ndlst, .fun=.getTxnyms)
  txnyms <- txnyms[1:length(txnyms)]
  names(txnyms) <- names(tree@ndlst)
  plyr::m_ply(1:length(tids), .pin)
  tree
}

# INPUT
data('mammals')
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
records.obj <- assembleRecords(records, taxonomy, tips=tree["tips"], .progress='time')
#load("records_obj.RData")
#save(records.obj, file="records_obj.RData")
lngs <- records.obj[['lineages']]
max_ages <- records.obj[['max.age']]
min_ages <- records.obj[['min.age']]
tids <- records.obj[['binomials']]

# recreate error
.pin <- function(i) {
  # unpack
  tid <- tids[i]
  print(tid)
  lng <- lngs[[i]]
  end <- ends[i]
  for(j in length(lng):1) {
    spns <- names(txnyms)[which(txnyms %in% lng[j])]
    spns <- spns[spns != tree@root]
    print(j)
    if(length(spns) == 0) {
      next
    }
    spns <- c(spns, unlist(sapply(spns, function(n) tree@ndlst[[n]][['ptid']])))
    rngs <- getSpnsAge(tree, ids=spns)
    bool <- rngs[ ,'start'] > end
    if(any(bool)) {
      rngs <- rngs[bool, ]
      rngs[rngs[ ,'end'] <= end, "end"] <- end
      # pinning is based on branch length
      prbs <- rngs$start - rngs$end
      e <- as.vector(sample(rngs$spn, prob=prbs, size=1))
      e_i <- which(rngs$spn == e)
      start <- runif(min=rngs$end[e_i], max=rngs$start[e_i], n=1)
      if(j != length(lng)) {
        tip_txnym <- lng[j+1]
      } else {
        tip_txnym <- lng[j]
      }
      pid <- paste0('p_', tid, sep='')
      tree <- addTip(tree, tid=tid, sid=e, start=start, end=end,
             pid=pid)
      tree@ndlst[[tid]][['txnym']] <- tip_txnym
      tree@ndlst[[pid]][['txnym']] <- lng[j]
      # add to txnyms list
      txnyms[[tid]] <<- tip_txnym
      # push out
      tree <<- tree
      break
    }
  }
}
.getTxnyms <- function(txnym, ...) {
  txnym
}
lngs <- records.obj[['lineages']]
max_ages <- records.obj[['max.age']]
min_ages <- records.obj[['min.age']]
tids <- records.obj[['binomials']]
ends <- sapply(1:length(min_ages),
               function(i) runif(min=min_ages[i], max=max_ages[i], n=1))

bool <- which(!tids %in% tree@tips)
tids <- tids[bool]
lngs <- lngs[bool]
ends <- ends[bool]
txnyms <- plyr::mlply(tree@ndlst, .fun=.getTxnyms)
txnyms <- txnyms[1:length(txnyms)]
names(txnyms) <- names(tree@ndlst)
plyr::m_ply(1:length(tids), .pin)

rndi <- which(tids %in% c("Nectomys_squamipes", "Kowalskia_sp."))

which(duplicated(tree@all))

spns <- getNdsSlt(tree, slt_nm='spn', tree@all)
bool <- spns < 0
which(bool)

kids <- getNdsKids(tree, ids=tree@nodes)
bool <- unlist(lapply(kids, function(k) any(duplicated(k))))
which(bool)

# PIN
outfolder <- file.path(output.dir, parent)

if (!file.exists(outfolder)) {
  dir.create(outfolder)
}
iterations <- 10
pin(tree, tids=tids, lngs=lngs, min_ages=min_ages,
    max_ages=max_ages, iterations=iterations, outfolder=outfolder)
