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
pin_is <- sort(as.numeric(sub('\\.RData', '', pinfiles)))
rndmfiles <- list.files(file.path(input_dir, paste0(parent, '_rndm')))
rnd_is <- sort(as.numeric(sub('\\.RData', '', rndmfiles)))

# CALCULATIONS IN PARALLEL
cat('Calculating ED by timeslice ....')
slice_dir <- file.path(output_dir, parent)
if(!file.exists(slice_dir)) {
  dir.create(slice_dir)
}
cat('.... real')
slice_files <- file.path(slice_dir, paste0('real_', pin_is, '.RData'))
is <- pin_is[!sapply(slice_files, file.exists)]
hldr <- foreach (i=is) %dopar% {
  tfl <- file.path(input_dir, paste0(parent, '_real'), paste0(i, '.RData'))
  sfl <- file.path(slice_dir, paste0('real_', i, '.RData'))
  cat ('\n........ [', i, ']', sep='')
  load(tfl)
  ed_slice <- calcEDBySlice(tree, time_cuts)
  save(ed_slice, file=sfl)
  rm(ed_slice)
}
cat('.... random')
slice_files <- file.path(slice_dir, paste0('rndm_', rnd_is, '.RData'))
is <- rnd_is[!sapply(slice_files, file.exists)]
hldr <- foreach (i=is) %dopar% {
  tfl <- file.path(input_dir, paste0(parent, '_rndm'), paste0(i, '.RData'))
  sfl <- file.path(slice_dir, paste0('rndm_', i, '.RData'))
  cat ('\n........ [', i, ']', sep='')
  load(tfl)
  ed_slice <- calcEDBySlice(tree, time_cuts)
  save(ed_slice, file=sfl)
  rm(ed_slice)
}
cat('Done.')

# TIMESTAMP
cat (paste0 ('\nSlicing finished at [', Sys.time (), ']'))