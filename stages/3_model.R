# D.J. Bennett
# 18/07/2015
# Model EDt=1 ~ EDt=0

# PARAMETERS
source('parameters.R')

# LIBS
library(ape)
library(nlme)

# DIRS
input_dir <- '2_slice'
output_dir <- '3_model'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

# INPUT
load(file=file.path(input_dir, 'ed_slices.RData'))
load(file=file.path(input_dir, 'nodedict.RData'))
load(file=file.path(input_dir, 'tree_code.RData'))
load(file=file.path(input_dir, 'constraint_trees.RData'))

# PROCESS
# loop through results from each tree
slopes <- res_slopes <- rep(NA, length(ed_slices))
for(i in 1:length(ed_slices)) {
  ed_slice <- log(ed_slices[[i]])
  # construct t0 and t1 by clade
  filler <- rep(NA, ncol(ed_slice))
  model_data <- data.frame(mean_ed=filler,
                            mean_diff=filler)
  rownames(model_data) <- colnames(ed_slice)
  for(j in 1:ncol(ed_slice)) {
    t0 <- ed_slice[-1,j]
    t1 <- ed_slice[-nrow(ed_slice),j]
    model_data[j, 'mean_ed'] <- mean(ed_slice[ ,j], na.rm=TRUE)
    model_data[j, 'mean_diff'] <- mean(t1-t0, na.rm=TRUE)
  }
  # clean data
  pull <- rowSums(model_data)
  pull <- !is.infinite(pull) & !is.na(pull)
  model_data <- model_data[pull, ]
  pull <- !constraint_trees[[i]]$tip.label %in% rownames(model_data)
  tree <- drop.tip(constraint_trees[[i]],
                   tip=constraint_trees[[i]]$tip.label[pull])
  model <- gls(mean_diff ~ mean_ed, data=model_data, method="ML",
                correlation=corPagel(value=1, phy=tree, fixed=TRUE))
  slopes[i] <- model$coefficients[2]
  model_data$residuals <- abs(model$residuals)
  model_res <- gls(residuals ~ mean_ed, data=model_data, method="ML",
                   correlation=corPagel(value=1, phy=tree, fixed=TRUE))
  res_slopes[i] <- model_res$coefficients[2]
}
# slope > 1: ED becomes more ED
# slope == 1: ED does not change
# slope == 0: ED is independent of past ED
# slope < 0: ED becomes less ED
rnd_lngs_slopes <- slopes[tree_code == 'rl']
rnd_ages_slopes <- slopes[tree_code == 'ra']
real_slopes <- slopes[tree_code == 'real']
sum(real_slopes > mean(rnd_ages_slopes, na.rm=TRUE)) / length(real_slopes)
sum(real_slopes > mean(rnd_lngs_slopes, na.rm=TRUE)) / length(real_slopes)
mean(real_slopes)
hist(real_slopes)
hist(rnd_ages_slopes)
hist(rnd_lngs_slopes)
# TODO -- calc heterosc
# res.slope > 0: high ED species act differently(i.e. oppositely to the model)
# res.slope < 0: low ED species act differently
# These results show the opposite of what I'd expect for a living fossil.
# They're more consistent with relicts, but perhaps phylogenetic scale will
# play a big role.
#save(slopes, res.slopes, model.data, file=file.path(output.dir, 'stats.Rd'))