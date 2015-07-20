# D.J. Bennett
# 18/07/2015
# Model EDt=1 ~ EDt=0

# DIRS
input.dir <- '2_slice'
output.dir <- '3_model'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
load (file=file.path (input.dir, 'ed_slices.Rd'))
load (file=file.path (input.dir, 'nodedict.Rd'))
ed.slices <- res
rm (res)
trees <- read.tree (file.path (input.dir, 'constraint_trees.tre'))

# PROCESS
# loop through results from each tree
slopes <- res.slopes <- rep (NA, length (ed.slices))
for (i in 1:length (ed.slices)) {
  ed.slice <- log (ed.slices[[i]])
  # construct t0 and t1 by clade
  filler <- rep (NA, ncol (ed.slice))
  model.data <- data.frame (mean.ed=filler,
                            mean.diff=filler)
  rownames (model.data) <- colnames (ed.slice)
  for (j in 1:ncol (ed.slice)) {
    t0 <- ed.slice[-1,j]
    t1 <- ed.slice[-nrow (ed.slice),j]
    model.data[j, 'mean.ed'] <- mean (ed.slice[ ,j], na.rm=TRUE)
    model.data[j, 'mean.diff'] <- mean (t1-t0, na.rm=TRUE)
  }
  # clean data
  pull <- rowSums (model.data)
  pull <- !is.infinite (pull) & !is.na (pull)
  model.data <- model.data[pull, ]
  pull <- !trees[[i]]$tip.label %in% rownames (model.data)
  tree <- drop.tip (trees[[i]], tip=trees[[i]]$tip.label[pull])
  model <- gls (mean.diff ~ mean.ed, data=model.data, method="ML",
                correlation=corPagel(value=1, phy=tree, fixed=TRUE))
  slopes[i] <- model$coefficients[2]
  model.data$residuals <- abs (model$residuals)
  model.res <- gls (residuals ~ mean.ed, data=model.data, method="ML",
                    correlation=corPagel(value=1, phy=tree, fixed=TRUE))
  res.slopes[i] <- model.res$coefficients[2]
}
# slope > 1: ED becomes more ED
# slope == 1: ED does not change
# slope == 0: ED is independent of past ED
# slope < 0: ED becomes less ED
sum (slopes > 0) / length (slopes)
mean (slopes)
hist (slopes)
t.test(slopes)
# TODO -- calc heterosc
# res.slope > 0: high ED species act differently (i.e. oppositely to the model)
# res.slope < 0: low ED species act differently
sum (res.slopes > 0) / length (slopes)
mean (res.slopes)
hist (res.slopes)
t.test (res.slopes)
# These results show the opposite of what I'd expect for a living fossil.
# They're more consistent with relicts, but perhaps phylogenetic scale will
# play a big role.