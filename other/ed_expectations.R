# FIGURE FOR SEEING HOW WE EXPECT ED t0 AND t1 WILL LOOK
x <- y <- NULL
newPlot <- function() {
  plot(x, x, type = 'n', xlab = '', ylab = '', xaxt = 'n',
       yaxt = 'n')
  lines(x, x, lty = 2, col = 'grey')
}
par(mfrow = c(2,2), mar = c(1, 0.5, 1, 0.5) + 0.1)
x <- 1:10
# panchronic
newPlot()
y <- c(0.5, 1, 1, 2, 5, 7, 9, 9, 11, 11)
lo <- loess(formula = y~x)
lines(predict(object = lo), col = 'red', lwd = 2)
# relict
newPlot()
y <- c(0, 3, 4, 6, 7, 8.5, 9.5, 9.5, 9.5, 9.9)
lo <- loess(formula = y~x)
lines(predict(object = lo), col = 'red', lwd = 2)
# relict + fuse
newPlot()
y <- c(1, 6, 7, 8, 9, 9.5, 9.5, 9.5, 6, 1)
lo <- loess(formula = y~x)
lines(predict(object = lo), col = 'red', lwd = 2)
# relict + panchronic
newPlot()
y <- c(5, 4, 3, 2, 1, 4, 6, 8, 10.5, 10.75)
lo <- loess(formula = y~x)
lines(predict(object = lo), col = 'red', lwd = 2)
