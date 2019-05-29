# Libs ----
library(foreach)
library(doMC)

# Functions ----
folder_gen <- iterate <- makeMdlData <- NULL
source(file.path('additional_analysis', 'bd_simulate_tools.R'))

# Vars ----
ncuts <- 5
niterations <- 10
ntips <- 1000
ncpus <- 4

# Dirs ----
outdir <- 'additional_analysis'
bd_dir <- folder_gen(file.path(outdir, 'bd'))
pan_dir <- folder_gen(file.path(outdir, 'pan'))
de_dir <- folder_gen(file.path(outdir, 'de'))
pf_dir <- folder_gen(file.path(outdir, 'pf'))

# Simulate ----
registerDoMC(ncpus)
cat('...bd\n')
tree_simulate(type = 'bd', flpth = bd_dir, overwrite = FALSE)
cat('...pan\n')
tree_simulate(type = 'pan', flpth = pan_dir, overwrite = FALSE)
cat('...de\n')
tree_simulate(type = 'de', flpth = de_dir, overwrite = FALSE)
cat('...pf\n')
tree_simulate(type = 'pf', flpth = pf_dir, overwrite = FALSE)

# Slice ----
slice(flpth = bd_dir)
slice(flpth = pan_dir)
slice(flpth = de_dir)
slice(flpth = pf_dir)

# Extract slices ----
bd_data <- makeMdlData(flpth = bd_dir)
pan_data <- makeMdlData(flpth = pan_dir)
de_data <- makeMdlData(flpth = de_dir)
pf_data <- makeMdlData(flpth = pf_dir)

# Save ----
save(bd_data, pan_data, de_data, pf_data,
     file = file.path('additional_analysis', 'sim_results.RData'))
