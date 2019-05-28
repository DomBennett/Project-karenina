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

# Plot ----
library(ggplot2)
all_data <- rbind(bd_data, pan_data, de_data, pf_data)
all_data$type <- c(rep('bd', nrow(bd_data)),
                   rep('pan', nrow(pan_data)),
                   rep('de', nrow(de_data)),
                   rep('pf', nrow(pf_data)))
all_data$id <- paste0(all_data$type, '_', all_data$id)
# drop multiples
all_data <- all_data[duplicated(all_data$id), ]
# individual
ggplot(all_data, aes(x = t0, y = t1)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1) +
  facet_grid(~type) + theme(legend.position = 'none')
# linear
ggplot(all_data, aes(x = t0, y = t1, colour = type)) +
  geom_smooth(method = 'loess') + geom_abline(slope = 1)
