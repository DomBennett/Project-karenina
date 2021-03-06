# D.J. Bennett
# 15/07/2015
# Pin fossils to phylogenies using PBDB

# TIMESTAMP
cat(paste0('\nPinning started at [', Sys.time(), ']'))

# LIBS
library(foreach)
library(doMC)
library(treeman)
pinParallelRand <- pinParallel <- assembleRecords <- NULL
source(file.path('tools', 'pin_tools.R'))

# PARAMETERS
overwrite <- parent <- treefile <- NULL
source('parameters.R')

# DIRS
input.dir <- '0_data'
output.dir <- '1_pin'
if (!file.exists(output.dir)) {
  dir.create(output.dir)
}

# READ TREE
cat('\nResearching tree ....')
tree <- readTree(file.path(input.dir, treefile), wndmtrx = FALSE)
tree_age <- getAge(tree)
tree_tips <- tree['tips']
# ----
# Reviewer comment: how many ranks?
# nms <- sample(tree_tips, 1000)
# res <- taxaResolve(nms, cache = TRUE, parent = parent)
# res <- res[!is.na(res[['score']]), ]
# ranks <- strsplit(x = res[['rank']], split = '\\|')
# mammal_ranks <- lapply(X = ranks, FUN = function(x) {
#   x[which(x == 'class'):length(x)]
# })
# length(unique(unlist(mammal_ranks))) - 2
# ----
txnyms <- searchTxnyms(tree, cache = TRUE, parent = parent)
wtxnyms <- sum(!is.na(txnyms))
# give parent ID to those with missing txnym
for (i in which(is.na(txnyms))) {
  nid <- names(txnyms)[i]
  while (is.na(txnyms[[i]])) {
    nid <- getNdSlt(tree, 'prid', nid)
    txnyms[[i]] <- txnyms[[nid]]
  }
}
tree <- setTxnyms(tree, txnyms)
cat('\nDone. [', wtxnyms, '/', tree['nall'], '] nodes with taxonyms.', sep = '')

# PALEODB
cat('\nRetrieving records ....')
fossilfile <- paste0(parent, '_records.csv')
if (overwrite | !file.exists(file.path(output.dir, fossilfile))) {
  # get
  cat('\n.... searching PBDB')
  records <- paleobioDB::pbdb_occurrences(limit = 'all', base_name = 'Aves',
                                          vocab = "pbdb",
                                          show = c("phylo", "ident"))
  # write out
  write.csv(records, file = file.path(output.dir, fossilfile),
             row.names = FALSE)
  # ensure factors are characters
  i <- sapply(records, is.factor)
  records[i] <- lapply(records[i], as.character)
} else {
  cat('\n.... pre-loading file')
  records <- read.csv(file = file.path(output.dir, fossilfile),
                      stringsAsFactors = FALSE)
}
cat('\nDone.')

# EXTRACT NAMES + LINEAGES
cat('\nAssembling records for pinning ....')
cat('.... [', nrow(records), '] records downloaded\n', sep = '')
# specify ranks in linages for pinning
taxonomy <- c('phylum', 'class', 'order', 'family', 'genus_name',
              'species_name')
records.obj <- assembleRecords(records, taxonomy, tree_tips = tree_tips,
                               .progress = 'time')
lineages <- records.obj[['lineages']]
max_age <- records.obj[['max.age']]
min_age <- records.obj[['min.age']]
binomials <- records.obj[['binomials']]
# drop fossils older than tree
pull <- max_age > tree_age
binomials <- binomials[!pull]
lineages <- lineages[!pull]
max_age <- max_age[!pull]
min_age <- min_age[!pull]
# drop records with too little lineage info
nlngs <- sapply(lineages, length)
pull <- nlngs > 3  # above family
binomials <- binomials[pull]
lineages <- lineages[pull]
max_age <- max_age[pull]
min_age <- min_age[pull]
# make sure all are within parent
pull <- sapply(lineages, function(x) ifelse(parent == 'Catarrhini',
                                            'Primates', parent) %in% x)
binomials <- binomials[pull]
lineages <- lineages[pull]
max_age <- max_age[pull]
min_age <- min_age[pull]
rm(records.obj)
rm(records)
cat('Done. Discovered [', length(binomials), '] records.', sep = '')

# FOSSIL STATS
cat('Determining fossil stats....\n')
fmls <- sapply(lineages, function(x) x[[3]])
nfmls <- length(unique(fmls))
cat('.... [', nfmls, '] families\n', sep = '')
gnra <- sapply(lineages, function(x) x[[length(x)]])
ngnra <- length(unique(gnra))
cat('.... [', ngnra, '] genera\n', sep = '')
mdpnts <- ((max_age - min_age)/2) + min_age
#hist(mdpnts, main = '', xlab = 'Fossil midpoint (MYA)')
cat('.... fossil midpoints:\t')
print(quantile(mdpnts))
cat('Done.\n')

# PIN
cat('\nPinning ....')
registerDoMC(ncpus)
cat('\n.... real')
pinfolder <- file.path(output.dir, paste0(parent, '_real'))
ntips <- pinParallel(tree, tids = binomials, lngs = lineages,
                     min_ages = min_age, max_ages = max_age,
                     pinfolder = pinfolder, tree_age = tree_age)
cat('\n.... random')
pinfolder <- file.path(output.dir, paste0(parent, '_rndm'))
ntips <- pinParallelRand(tree, tids = binomials, pinfolder = pinfolder,
                         min_ages = min_age, tree_age = tree_age)
cat('\nDone.')

# TIMESTAMP
cat(paste0('\nPinning finished at [', Sys.time(), ']'))
