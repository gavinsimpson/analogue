##' @title Number of samples per gradient segments
##'
##' @description The number of samples in sections along the gradient
##' is a useful diagnostic as to the quality of reconstructions at
##' gradient values within those sections.
##'
##' @details The sampling design of a training set, i.e. the number of
##' samples taken at points along the gradient, can influence the
##' uncertainty in the transfer function predictions at those values
##' of the gradient. Poorly sampled sections of the gradient may have
##' far larger RMSEP than the overall model RMSEP.
##'
##' @param grad numeric; vector of gradient values
##' @param n numeric; number of segments to partition the gradient into
##'
##' @return Numeric vector of length \code{n} containing the numbers of
##' samples per gradient segment.
##'
##' @author Gavin L. Simpson
##'
##' @keywords utilities
##'
##' @examples
##' data(SumSST)
##' evenSample(SumSST)  ## not an even sample...
`evenSample` <- function(grad, n = 10) {
    segs <- cut(grad, breaks = n) ##, labels = paste0("seg", seq_len(n)))
    Nseg <- tapply(grad, segs, length)
    attr(Nseg, "gradient") <- deparse(substitute(grad))
    attr(Nseg, "numSegments") <- n
    class(Nseg) <- "evenSample"
    Nseg
}
