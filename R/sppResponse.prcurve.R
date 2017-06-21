`sppResponse.prcurve` <- function(x, n = 100, ...) {
    Pred <- function(x, n) {
        rn <- range(x$lambda)
        newx <- seq(from = rn[1], to = rn[2], length = n)
        isGAM <- inherits(x$model, "gam")
        if (isGAM) {
            newx <- data.frame(lambda = newx)
            fitted.values <- predict(x$model, newdata = newx, type = "response")
            fitted.values <- cbind(newx, as.vector(fitted.values))
        } else {
            fitted.values <- predict(x$model, newx)
        }
    observed <- x[c("lambda","x")]
    names(fitted.values) <- names(observed) <- c("gradient","response")
    list(observed = observed, fitted.values = fitted.values)
  }
  out <- lapply(x$smooths, Pred, n = n)
  class(out) <- "sppResponse"
  out
}
