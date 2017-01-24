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
source(file.path('tools', 'slice_tools_tm.R'))

# LIST FILES
pinfiles <- list.files(file.path(input_dir, paste0(parent, '_real')))
rndlngfiles <- list.files(file.path(input_dir, paste0(parent, '_rndm_lngs')))
rndagefiles <- list.files(file.path(input_dir, paste0(parent, '_rndm_ages')))
tree_files <- c(file.path(input_dir, paste0(parent, '_real'), pinfiles),
                file.path(input_dir, paste0(parent, '_rndm_lngs'), rndlngfiles),
                file.path(input_dir, paste0(parent, '_rndm_ages'), rndagefiles))
tree_code <- c(rep('p', each=length(pinfiles)),
               rep('rl', each=length(rndlngfiles)),
               rep('ra', each=length(rndagefiles)))
rm(pinfiles, rndagefiles, rndlngfiles)

# CALCULATIONS IN PARALLEL
cat('Calculating ED by timeslice ....')
ed_slices <- foreach (i=1:length(tree_files)) %dopar% {
  cat ('\n........ [', i, ']', sep='')
  load(tree_files[[i]])
  calcEDBySlice(tree, time_cuts)
}
cat('Done.')

# OUTPUT
save(tree_code, file=file.path(output_dir, 'tree_code.RData'))
save(ed_slices, file=file.path(output_dir, 'ed_slices.RData'))

# TIMESTAMP
cat (paste0 ('\nSlicing finished at [', Sys.time (), ']'))