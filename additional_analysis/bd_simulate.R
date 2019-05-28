# Libs ----
library(foreach)
library(doMC)

# Functions ----
folder_gen <- iterate <- makeMdlData <- NULL
source(file.path('additional_analysis', 'bd_simulate_tools.R'))

# Vars ----
ncuts <- 10
niterations <- 2
ntips <- 1000
ncpus <- 2

# Dirs ----
outdir <- 'additional_analysis'
bd_dir <- folder_gen(file.path(outdir, 'bd'))
pan_dir <- folder_gen(file.path(outdir, 'pan'))
eph_dir <- folder_gen(file.path(outdir, 'eph'))
de_dir <- folder_gen(file.path(outdir, 'de'))
pf_dir <- folder_gen(file.path(outdir, 'pf'))

# Simulate ----
registerDoMC(ncpus)
cat('...bd\n')
iterate(type = 'bd', flpth = bd_dir, overwrite = FALSE)
cat('...pan\n')
iterate(type = 'pan', flpth = pan_dir, overwrite = FALSE)
cat('...eph\n')
iterate(type = 'eph', flpth = eph_dir, overwrite = FALSE)
cat('...de\n')
iterate(type = 'de', flpth = de_dir, overwrite = FALSE)
cat('...pf\n')
iterate(type = 'pf', flpth = pf_dir, overwrite = FALSE)

# Extract slices ----
# bd
ed_files <- list.files(path = bd_dir, pattern = '.RData')
bd_data <- makeMdlData(ed_files = file.path(bd_dir, ed_files))
# pan
ed_files <- list.files(path = pan_dir, pattern = '.RData')
pan_data <- makeMdlData(ed_files = file.path(pan_dir, ed_files))
# eph
ed_files <- list.files(path = eph_dir, pattern = '.RData')
eph_data <- makeMdlData(ed_files = file.path(eph_dir, ed_files))
# de
ed_files <- list.files(path = de_dir, pattern = '.RData')
de_data <- makeMdlData(ed_files = file.path(de_dir, ed_files))
# pf
ed_files <- list.files(path = pf_dir, pattern = '.RData')
pf_data <- makeMdlData(ed_files = file.path(pf_dir, ed_files))

# Plot ----
library(ggplot2)
all_data <- rbind(bd_data, pan_data, eph_data, de_data, pf_data)
all_data$type <- c(rep('bd', nrow(bd_data)),
                   rep('pan', nrow(pan_data)),
                   rep('eph', nrow(eph_data)),
                   rep('de', nrow(de_data)),
                   rep('pf', nrow(pf_data)))
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
  geom_smooth()

