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
registerDoMC(ncpus)
library(treeman)
source(file.path('tools', 'slice_tools.R'))

# LIST FILES
pinfiles <- list.files(file.path(input_dir, paste0(parent, '_real')))
rndmfiles <- list.files(file.path(input_dir, paste0(parent, '_rndm')))
tree_files <- c(file.path(input_dir, paste0(parent, '_real'), pinfiles),
                file.path(input_dir, paste0(parent, '_rndm'), rndmfiles))
tree_code <- c(rep('pin', each=length(pinfiles)),
               rep('rnd', each=length(rndmfiles)))
rm(pinfiles, rndmfiles)

# CALCULATIONS IN PARALLEL
cat('Calculating ED by timeslice ....')
slice_dir <- file.path(output_dir, parent)
if(!file.exists(slice_dir)) {
  dir.create(slice_dir)
}
slice_files <- file.path(slice_dir,paste0(1:length(tree_files), '.RData'))
is <- which(!sapply(slice_files, file.exists))
foreach (i=is) %dopar% {
  cat ('\n........ [', i, ']', sep='')
  load(tree_files[[i]])
  ed_slice <- calcEDBySlice(tree, time_cuts)
  save(ed_slice, file=slice_files[i])
  rm(ed_slice)
}
cat('Done.')

# OUTPUT
save(tree_code, file=file.path(slice_dir, 'tree_code.RData'))

# TIMESTAMP
cat (paste0 ('\nSlicing finished at [', Sys.time (), ']'))