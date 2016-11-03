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
  rand.slice <- ed.slices[[sample (1:100, 1)]]
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