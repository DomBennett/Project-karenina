# D.J. Bennett
# 18/07/2015
# Take time slices of time-callibrated trees

# PARAMETERS
time.cuts <- 10

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
trees <- read.tree (file.path (input.dir, 'pinned.tre'))

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
save (res, file=file.path (output.dir, 'ed_slices.Rd'))
save (nodedict, file=file.path (output.dir, 'nodedict.Rd'))
write.tree (constraint.trees,
            file=file.path (output.dir, 'constraint_trees.tre'))