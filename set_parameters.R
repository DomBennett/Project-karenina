# D.J. Bennett
# 27/07/2015
# Set parameters

treefile <- 'bininda_mammalia.tre'  # name of tree file in 0_data/
parent <- "catarrhini"  # name of parent clade for species in tree
iterations <- 100  # number of pinned trees
overwrite <- FALSE  # overwrite download files
time.cuts <- c (1.29985, 3.9605, 14.1815, 28.465, 44.95,
                61, 83.25, 122.75, 154.25, 168.8, 187.7,
                219.15, 242.1, 249.7, 256, 266.05, 285.6)
# ages to cut tree (here using middle epoch values as
# these are the units of sampling in the pdbd)
ncpus <- 2  # make 1 if no parallisation (UNIX-only)
save (treefile, parent, iterations, overwrite, time.cuts,
      ncpus, file='parameters.Rd')

###########################
# Middle values of epochs #
###########################
# Pleistocene (0.0117 - 2.588, 1.29985)
# Pliocene (2.588 - 5.333, 3.9605)
# Miocene (5.333 - 23.03, 14.1815)
# Oligocene (23.03 - 33.9, 28.465)
# Eocene (33.9 - 56.0, 44.95)
# Paleocene (56.0 - 66.0, 61)
# Cretaceous Upper (66.0 - 100.5, 83.25)
# Cretaceous Lower (100.5 - 145.0, 122.75)
# Jurrasic Upper (145.0 - 163.5, 154.25)
# Jurrasic Middle (163.5 - 174.1, 168.8)
# Jurrasic Lower (174.1 - 201.3, 187.7)
# Triassic Upper (201.3 - 237.0, 219.15)
# Triassic Middle (237.0 - 247.2, 242.1)
# Triassic Lower (247.2 - 252.2, 249.7)
# Permian Lopingian (252.2 - 259.8, 256)
# Permian Guadaulupian (259.8 - 272.3, 266.05)
# Permian Cisuralian (272.3 - 298.9, 285.6)