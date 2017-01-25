# FUNCTIONS
slice <- function(tree, tm_ct) {
  getTx <- function(nd) {
    txnyms <- nd[['txnym']]
    txnyms[length(txnyms)]
  }
  tree_age <- getAge(tree)
  nd_spns <- getSpnsAge(tree, tree@all, tree_age=tree_age)
  nd_spns[['spn']] <- as.character(nd_spns[['spn']])
  # drop tips extinct by time cut
  tp_ages <- nd_spns[nd_spns[['spn']] %in% tree@tips, c('spn', 'end')]
  bool <- tp_ages[['end']] > tm_ct
  if(any(bool)) {
    to_drp <- tp_ages[['spn']][bool]
    tree <- rmTips(tree, tids=to_drp)
    nd_spns <- nd_spns[nd_spns[['spn']] %in% tree@all, ]
  }
  # slice tree at interval
  slcd <- timeslice(tree=tree, tm_ct=tm_ct,
                    nd_spns=nd_spns)
  # replace IDs with txnyms
  txnyms <- sapply(slcd@ndlst[slcd@tips], getTx)
  new_ids <- slcd@tips
  bool <- grepl('^n[0-9]', new_ids)
  new_ids[bool] <- txnyms[bool]
  duplicates <- which(duplicated(new_ids))
  if(length(duplicates) > 0) {
    for(i in duplicates) {
      bool <- new_ids[i] == new_ids
      if(sum(bool) == 1) {
        next
      }
      ncps <- sum(bool)
      new_ids[bool] <- paste0(new_ids[bool], '_', 1:ncps)
    }
  }
  setNdsID(slcd, slcd@tips, new_ids)
}

# PARAMETERS
parent <- 'hominoidea'

# DIRS
input_dir <- '1_pin'
output_dir <- '2_slice'
if(!file.exists(output_dir)) {
  dir.create(output_dir)
}

# LIBS
library(treeman)
source(file.path('tools', 'slice_tools_tm.R'))

# LIST FILES
pinfile <- list.files(file.path(input_dir, paste0(parent, '_real')))[1]
tree_file <- file.path(input_dir, paste0(parent, '_real'), pinfile)

# LOAD
load(tree_file)

# SLICE
# Miocene
slcd_1 <- slice(tree, 14.1815)
# pleisotcene
slcd_2 <- slice(tree, 1.29985)

# PLOT
pdf(file='reconstructions.pdf', w=7, h=15)
par(mfrow=c(3,1))
plot(as(tree, 'phylo'), main='Eocene - Present', cex=0.5)
ape::axisPhylo()
par(mar=c(2, 0.1, 1, 0.1))
plot(as(slcd_1, 'phylo'), main='Miocene (14 MYA)', cex=1.1)
plot(as(slcd_2, 'phylo'), main='Pleistocene (1.3 MYA)', cex=1.1)
dev.off()