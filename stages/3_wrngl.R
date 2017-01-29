# D.J. Bennett
# 29/01/2017

# TIMESTAMP
cat (paste0 ('\nWrangling started at [', Sys.time (), ']'))

# PARAMETERS
source('parameters.R')

# LIBS
source(file.path ('tools', 'wrngl_tools.R'))

# DIRS
input_dir <- '2_slice'
output_dir <- '3_wrngl'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_dir <- file.path(input_dir, parent)

# INPUT
load(file=file.path(input_dir, 'tree_code.RData'))

# MAKE MODEL DATA
cat('Merging slices in model data ....\n')
ed_files <- list.files(input_dir)
ed_files <- ed_files[ed_files != 'tree_code.RData']
ed_files <- file.path(input_dir, ed_files)
mdl_data <- data.frame(t0=NA, t1=NA, tmsplt=NA, id=NA, cnt=NA)
for(i in 1:length(ed_files)) {
  cat('....[', i, ']\n', sep='')
  load(ed_files[[i]])
  tmp <- sapply(1:(nrow(ed_slice) - 1), extrct)
  rm(ed_slice)
}
mdl_data <- mdl_data[-1, ]
cat('Done.\n')

# CHECK
cat('Checking....\n')
tmsplts <- unique(mdl_data$tmsplt)
chck <- rep(NA, length(tmsplts))
for(i in 1:length(tmsplts)) {
  chck[i] <- sum(duplicated(mdl_data$id[mdl_data$tmsplt == tmsplts[i]]))
}
pdf(file.path(output_dir, paste0(parent, '_check.pdf')))
hist(mdl_data$cnt)
plot(t0~cnt, data=mdl_data)
dev.off()
if(all(chck == 0)) {
  cat('Done. Check [passed].\n')
} else {
  cat('Done. Check [failed].\n')
}

# SORT
cat('Sorting data....\n')
tmsplts <- c("1.29985-0", "3.9605-1.29985", "14.1815-3.9605",
             "28.465-14.1815", "44.95-28.465", "61-44.95",
             "83.25-61", "122.75-83.25", "154.25-122.75")
epochs <- c('Pe-Re', 'Pi-Pe', 'Mi-Pi', 'Ol-Mi', 'Eo-Ol',
            'Pa-Eo', 'CU-Pa', 'CL-CU', 'JU-CL')
mtchng <- match(mdl_data$tmsplt, tmsplts)
mdl_data$epoch <- epochs[mtchng]
mdl_data$epoch <- factor(mdl_data$epoch,
                          levels=epochs)
# add tree info
library(treeman)
tree <- readTree(file.path('0_data', treefile),
                 wndmtrx=FALSE)
# identify all fossil record nds
mdl_data$fssl_nd <- !mdl_data$id %in% tree['all']
mdl_data$fssl_nd <- factor(mdl_data$fssl_nd)
rm(tree)
# TODO: add taxonomic using pinned trees
cat('Done.\n')

# SAVE
save(mdl_data, file=file.path(output_dir, paste0(parent, '.RData')))

# TIMESTAMP
cat (paste0 ('\nWrangling finished at [', Sys.time (), ']'))