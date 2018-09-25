# D.J. Bennett
# 29/01/2017

# TIMESTAMP
cat(paste0('\nWrangling started at [', Sys.time(), ']'))

# LIBS
library(treeman)
getOrdrAndGenera <- makeMdlData <- NULL
source(file.path('tools', 'wrngl_tools.R'))

# PARAMETERS
treefile <- parent <- NULL
source('parameters.R')

# DIRS
input_dir <- '2_slice'
output_dir <- '3_wrngl'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_dir <- file.path(input_dir, parent)

# MAKE MODEL DATA
cat('Merging slices in model data ....\n')
fls <- list.files(input_dir, pattern = '.RData')
# real
cat('\n.... real')
ed_files <- fls[grepl('^real_', fls)]
ed_files <- file.path(input_dir, ed_files)
mdl_data <- makeMdlData(ed_files)
# random
cat('\n.... random')
ed_files <- fls[grepl('^rndm_', fls)]
ed_files <- file.path(input_dir, ed_files)
rnd_data <- makeMdlData(ed_files)
cat('Done.\n')

# CHECK
cat('Checking....\n')
tmsplts <- unique(mdl_data$tmsplt)
chck <- rep(NA, length(tmsplts))
for (i in 1:length(tmsplts)) {
  chck[i] <- sum(duplicated(mdl_data$id[mdl_data$tmsplt == tmsplts[i]]))
}
pdf(file.path(output_dir, paste0(parent, '_check.pdf')))
hist(mdl_data$cnt)
plot(t0~cnt, data = mdl_data)
dev.off()
if (all(chck == 0)) {
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
# real
mtchng <- match(mdl_data$tmsplt, tmsplts)
mdl_data$epoch <- epochs[mtchng]
mdl_data$epoch <- factor(mdl_data$epoch, levels = epochs)
# rand
mtchng <- match(rnd_data$tmsplt, tmsplts)
rnd_data$epoch <- epochs[mtchng]
rnd_data$epoch <- factor(rnd_data$epoch, levels = epochs)
# add tree info
tree <- readTree(file.path('0_data', treefile), wndmtrx = FALSE)
# identify all fossil record nds
mdl_data$fssl_nd <- !mdl_data$id %in% tree['all']
mdl_data$fssl_nd <- factor(mdl_data$fssl_nd)
rm(tree)
flp <- file.path('1_pin', paste0(parent, '_real'))
load(file.path(flp, sample(list.files(flp, pattern = '.RData'), 1)))
# extract genus and order(ish)
orders <- c('Monotremata', 'Proboscidea', 'Sirenia', 'Tenrecidae',
            'Macroscelidea', 'Hyracoidea', 'Chrysochloridae', 'Scandentia',
            'Primates', 'Dermoptera', 'Rodentia', 'Lagomorpha', 'Carnivora',
            'Cetartiodactyla', 'Chiroptera', 'Insectivora', 'Perissodactyla',
            'Pholidota', 'Litopterna', 'Xenarthra', 'Notoungulata',
            'Metatheria')
tips <- tree['tips']
genera <- unique(sub('_.*$', '', tips))
res <- plyr::mdply(tree['all'], .fun = getOrdrAndGenera, .progress = 'time',
                   tree = tree, orders = orders, genera)
# add to data
mtchng <- match(mdl_data[['id']], res[['nid']])
mdl_data[['order']] <- res[mtchng, 'order']
mdl_data[['genus']] <- res[mtchng, 'genus']
cat('Done.\n')

# SAVE
save(mdl_data, rnd_data, file = file.path(output_dir, paste0(parent, '.RData')))

# TIMESTAMP
cat(paste0('\nWrangling finished at [', Sys.time(), ']'))
