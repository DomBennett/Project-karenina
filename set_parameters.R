# D.J. Bennett
# 27/07/2015
# Set parameters

treefile <- 'bininda_mammalia.tre'  # name of tree file in 0_data/
parent <- "catarrhini"  # name of parent clade for species in tree
iterations <- 100  # number of pinned trees
overwrite <- FALSE  # overwrite download files
time.cuts <- 10  # number of cuts in tree
ncpus <- 2  # make 1 if no parallisation (UNIX-only)
save (treefile, parent, iterations, overwrite, time.cuts,
      ncpus, file='parameters.Rd')
