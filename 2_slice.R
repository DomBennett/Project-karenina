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
nodedict <- list ()
for (i in 1:length (trees)) {
  part.res <- calcEDBySlice (trees[[i]], time.cuts)
  res <- c (res, list (part.res))
  nnode <- getSize (trees[[i]]) + trees[[i]]$Nnode
  for (j in 1:nnode) {
    node.name <- paste0 ('n', j)
    nodedict[node.name] <- list (getChildren (trees[[i]], j))
  }
}

# OUTPUT
save (res, file=file.path (output.dir, 'ed_slices.Rd'))
save (nodedict, file=file.path (output.dir, 'nodedict.Rd'))