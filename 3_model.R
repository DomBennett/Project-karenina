# D.J. Bennett
# 18/07/2015
# Model EDt=1 ~ EDt=0
# TODO: look up how to build correlation matrix

# DIRS
input.dir <- '2_slice'
output.dir <- '3_model'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
load (file=file.path (input.dir, 'ed_slices.Rd'))
#load (file=file.path (input.dir, 'nodedict.Rd'))
ed.slices <- res
rm (res)

# PROCESS
# loop through results from each tree
slopes <- rep (NA, length (ed.slices))
for (i in 1:length (ed.slices)) {
  ed.slice <- ed.slices[[i]]
  # construct t0 and t1 by clade
  t1 <- t0 <- NULL
  for (j in 1:ncol (ed.slice)) {
    t1 <- c (t1, ed.slice[-1,j])
    t0 <- c (t0, ed.slice[-nrow (ed.slice),j])
  }
  #plot (t1 ~ t0,
  #      xlab=expression('ED'['t=0']),
  #      ylab=expression('ED'['t=1']),)
  #abline (a=0, b=1)
  model.data <- data.frame (t1, t0)
  pull <- rowSums (model.data)
  pull <- !is.infinite (pull) & !is.na (pull)
  model.data <- model.data[pull, ]
  #cor.p <- corPagel(1, pin.res, fixed=FALSE)
  model <- gls (t1 ~ t0, data=model.data, method="ML")
  slopes[i] <- model$coefficients[2]
}
sum (slopes > 1)/length (slopes)
mean (slopes)
