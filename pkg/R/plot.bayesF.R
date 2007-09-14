## Plots probability of Analogue based on
## Bayes factors
plot.bayesF <- function(x,
                        xlab = NULL,
                        ylab = "Pr (A+ | d)",
                        col = "red",
                        abline.col = "lightgrey",
                        ...) {
  if(!inherits(x, "bayesF"))
    stop("Plot method only for objects of class \"bayesF\".")
  if(is.null(xlab))
     xlab <- paste("Dissimilarity (", x$method, ")", sep = "")
  prob.pos <- x$posterior$pos / (1 + x$posterior$pos)
  prob.pos[is.nan(prob.pos)] <- 1
  plot(rev(x$roc.points), prob.pos, type = "n",
       ylab = ylab, xlab = xlab, axes = FALSE)
  abline(v = x$optimal, lty = "dotted", col = abline.col)
  lines(rev(x$roc.points), prob.pos, col = col)
  axis(1)
  axis(2)
  box()
  invisible(x)
}
