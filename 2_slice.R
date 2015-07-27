# D.J. Bennett
# 18/07/2015
# Take time slices of time-callibrated trees

# TIMESTAMP
cat (paste0 ('\nSlicing started at [', Sys.time (), ']'))

# PARAMETERS
load ('parameters.Rd')

# DIRS
input.dir <- '1_pin'
output.dir <- '2_slice'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# LIBS
library (MoreTreeTools)
source (file.path ('tools', 'palaeo_tools.R'))

# INPUT
pinfile <- paste0 (parent, '_pinned.tre')
pin.trees <- read.tree (file.path (input.dir, pinfile))
randfile <- paste0 (parent, '_rand.tre')
rand.trees <- read.tree (file.path (input.dir, randfile))
trees <- c (pin.trees, rand.trees)
rand.bool <- c (rep (FALSE, length (pin.trees)),
                rep (TRUE, length (rand.trees)))

# CALC ED BY TIMESLICE
res <- list ()  # list of matrices
nodedicts <- list ()
constraint.trees <- list ()
for (i in 1:length (trees)) {
  # get ED by timeslice
  part.res <- calcEDBySlice (trees[[i]], time.cuts)
  res <- c (res, list (part.res))
  # get nodedict for tree
  nodedict <- getNodedict (trees[[i]])
  nodedicts <- c (nodedicts, list (nodedict))
  # get constraint tree
  constraint.tree <- addTips2Nodes (trees[[i]])
  constraint.trees <- c (constraint.trees,
                         list (constraint.tree))
}
class (constraint.trees) <- 'multiPhylo'

# OUTPUT
save (rand.bool, file=file.path (output.dir, 'randbool.Rd'))
save (res, file=file.path (output.dir, 'ed_slices.Rd'))
save (nodedict, file=file.path (output.dir, 'nodedict.Rd'))
write.tree (constraint.trees,
            file=file.path (output.dir, 'constraint_trees.tre'))

# TIMESTAMP
cat (paste0 ('\nSlicing finished at [', Sys.time (), ']'))