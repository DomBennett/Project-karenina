# D.J. Bennett
# 18/07/2015
# Take time slices of time-callibrated trees

# TIMESTAMP
cat (paste0 ('\nSlicing started at [', Sys.time (), ']'))

# PARAMETERS
source('parameters.R')

# DIRS
input_dir <- '1_pin'
output_dir <- '2_slice'
if(!file.exists(output_dir)) {
  dir.create(output_dir)
}

# LIBS
require (foreach)
require (doMC)
registerDoMC (ncpus)
library(ape)
source(file.path('tools', 'slice_tools.R'))

# INPUT
pinfile <- paste0(parent, '_real.tre')
pin_trees <- read.tree(file.path(input_dir, pinfile))
randfile <- paste0(parent, '_rndm_lngs.tre')
rndm_lngs_trees <- read.tree(file.path(input_dir, randfile))
randfile <- paste0(parent, '_rndm_ages.tre')
rndm_ages_trees <- read.tree(file.path(input_dir, randfile))
trees <- c (pin_trees, rndm_lngs_trees, rndm_ages_trees)
tree_code <- rep(c('real', 'rl', 'ra'), each=iterations)

# CALCULATIONS IN PARALLEL
cat('Calculating ED by timeslice ....')
ed_slices <- foreach (i=1:length(trees)) %dopar% {
  cat ('\n........ [', i, ']', sep='')
  calcEDBySlice(trees[[i]], time.cuts)
}
cat('Done.')
cat('Calculating nodedicts ....')
nodedict <- foreach (i=1:length(trees)) %dopar% {
  cat ('\n........ [', i, ']', sep='')
  # get nodedict for tree
  getNodedict(trees[[i]])
}
cat('Done.')
cat('Calculating constraint trees ....')
constraint_trees <- foreach (i=1:length(trees)) %dopar% {
  cat ('\n........ [', i, ']', sep='')
  # get constraint tree
  addTips2Nodes(trees[[i]])
}
class(constraint_trees) <- 'multiPhylo'
cat('Done.')

# OUTPUT
save(tree_code, file=file.path(output_dir, 'tree_code.RData'))
save(ed_slices, file=file.path(output_dir, 'ed_slices.RData'))
save(nodedict, file=file.path(output_dir, 'nodedict.RData'))
save(constraint_trees, file=file.path(output_dir, 'constraint_trees.RData'))

# TIMESTAMP
cat (paste0 ('\nSlicing finished at [', Sys.time (), ']'))