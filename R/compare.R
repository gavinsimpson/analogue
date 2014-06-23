##' @title Compare proxies across two data sets
##'
##' @description.. content for \description{} (no empty lines) ..
##'
##' @details To Do
##'
##' @param x data frame; training set samples to compare against
##' @param ... other arguments passed to methods.
##'
##' @return To Do
##'
##' @author Gavin L. Simpson
##'
##' @rdname compare
##'
 `compare` <- function(x, ...) {
    UseMethod("compare")
}

##' @param y data frame; passive or core samples
##' @param env numeric vector of environmental or contraint data for residual length ordination
##' @param ordination character; which constrained ordination method to use
##' @param method character; which dissimilarity method to use. See \code{distance}.
##' @param transform character: should a transformation be applied to the data. Ignored.
##'
##' @rdname compare
##'
`compare.default` <- function(x, y, env,
                              ordination = "rda",
                              method = "chord",
                              transform = NULL,
                              ...) {
    joint <- join(x, y, split = TRUE)
    x2 <- joint$x
    y2 <- joint$y
    n2.x <- n2(x2, which = "species")
    out <- data.frame(sumMissing  = rowSums(y2[, is.infinite(n2.x)]),
                      sumPoorOpt  = rowSums(y2[, n2.x <= 5]),
                      closestSamp = minDC(ana.fit <- analog(x2, y2))$minDC,
                      residLen    = residLen(x, env, y, method = ordination)[["passive"]])
    out
}
