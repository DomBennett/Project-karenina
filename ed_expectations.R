# FIGURE FOR SEEING HOW WE EXPECT ED t0 AND t1 WILL LOOK

newPlot <- function() {
  plot(x, x, type = 'n', xlab='', ylab='', xaxt='n',
       yaxt='n')
  lines(x, x, lty = 2, col='grey')
}
par(mfrow=c(2,2), mar=c(1, 0.5, 1, 0.5) + 0.1)
x <- 1:10
# panchronic
newPlot()
y <- c(1, 1, 1, 1, 2, 6, 10, 10.25, 10.5, 10.75)
lines(x, y)
# relict
newPlot()
y <- c(1, 6, 7, 8, 9, 9.5, 9.5, 9.5, 9.5, 9.5)
lines(x, y)
# relict + fuse
newPlot()
y <- c(1, 6, 7, 8, 9, 10, 10, 10, 6, 1)
lines(x, y)
# relict + panchronic
newPlot()
y <- c(10, 8, 6, 4, 1, 4, 6, 8, 10.5, 10.75)
lines(x, y)