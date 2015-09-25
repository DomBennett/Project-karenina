# D.J. Bennett
# 18/07/2015
# Model EDt=1 ~ EDt=0

# PARAMETERS
load ('parameters.Rd')

# LIBS
library (ape)
library (nlme)

# DIRS
input.dir <- '2_slice'
output.dir <- '3_model'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
load (file=file.path (input.dir, 'ed_slices.Rd'))
load (file=file.path (input.dir, 'nodedict.Rd'))
load (file=file.path (input.dir, 'randbool.Rd'))
ed.slices <- res
rm (res)
trees <- read.tree (file.path (input.dir, 'constraint_trees.tre'))

# New idea
# randomly select N clades from ED slice across all trees from one time slice
# N is calculated as the number of species in the smallest time slice
# Rank the N clades by their ED value
# Repeat for all time slices
# Repeat I times to generate a distirbution of estimates per time slice
# Model ED rank by time slice
# Use ED rank because this is comparable between time slices, straight ED would still
# be affected by sampling by time slice
# Question: do I still need to use a GLS? N may be too small.

y <- x <- NULL
for (j in 1:nrow (ed.slices[[1]])) {
  ty <- ed.slices[[1]][j, ]
  tx <- as.numeric (rep (rownames (ed.slices[[1]])[j], length (ty)))
  y <- c (y, ty)
  x <- c (x, tx)
}
pull <- is.na (y)
y <- y[!pull]
x <- x[!pull]

sample.n <- 3
t0 <- t1 <- NULL
for (i in 1:iterations) {
  # select random ed.slice
  rand.slice <- ed.slices[[sample (101:200, 1)]]
  # select random time points to compare
  rand.t0 <- sample (1:(nrow (rand.slice) - 1), 1)
  rand.t1 <- rand.t0 + 1
  # select sample.n lineages between t0 and t1 to compare
  possibles <- which (!is.na (rand.slice[rand.t0, ]) &
                        !is.na (rand.slice[rand.t1, ]))
  if (length (possibles) < sample.n) {
    next
  }
  rand.n <- sample (possibles, sample.n)
  # get their ranks
  t0.rank <- rank (rand.slice[rand.t0, rand.n])
  t1.rank <- rank (rand.slice[rand.t1, rand.n])
  # add to global list
  t0 <- c (t0, t0.rank)
  t1 <- c (t1, t1.rank)
}
plot (t1 ~ t0)
model <- lm (t1 ~ t0)
summary (model)
abline (model)

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
rand.slopes <- slopes[rand.bool]
pin.slopes <- slopes[!rand.bool]
sum (pin.slopes > mean (rand.slopes, na.rm=TRUE)) / length (pin.slopes)
mean (pin.slopes)
hist (pin.slopes)
hist (rand.slopes)
t.test (pin.slopes, rand.slopes)
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
save (slopes, res.slopes, model.data, file=file.path (output.dir, 'stats.Rd'))