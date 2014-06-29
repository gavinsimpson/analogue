###########################################################################
##                                                                       ##
## Internal functions for package analogue - not meant to be used by     ##
## users.                                                                ##
##                                                                       ##
###########################################################################

##' @title Cumulative weighted mean for a vector of distances
##'
##' @param weights numeric vector of weights
##' @param y numeric vector to calculate the weighted mean of
##' @param drop logical; drop spurious zero distance
##' @param kmax numeric; upper limit on number of analogues to include
##'
##' @return a numeric vector of length \code{kmax}.
##'
##' @author Gavin L. Simpson
cumWmean <- function(weights, y, drop = TRUE, kmax) {
    ## as weights are the distances, I could probably combine
    ## mean and weighted mean versions of this function?
    if(missing(kmax))
        kmax <- length(y)
    ##if (length(weights) != length(y))
    ##  stop("'y' and 'weights' must have the same length")
    nas <- is.na(weights)
    ord <- order(weights[!nas])
    if(drop) {
        weights <- 1 / weights[!nas][ord][-1]
        env <- y[!nas][ord][-1]
    } else {
        weights <- 1 / weights[!nas][ord]
        env <- y[!nas][ord]
    }
    K <- seq_len(kmax)
    cumsum(weights[K] * env[K]) / cumsum(weights[K])
}

##' @title Cumulative mean for a vector of distances
##'
##' @param dis the distances to sort by
##' @param y the vector of values to calculate mean of
##' @param drop logical; drop spurious zero distance
##' @param kmax numeric; upper limit on number of analogues to include
##'
##' @return a numeric vector of length \code{kmax}.
##'
##' @author Gavin L. Simpson
cummean <- function(dis, y, drop = TRUE, kmax) {
    if(missing(kmax))
        kmax <- length(y)
    nas <- is.na(dis)
    ord <- order(dis[!nas])
    y <- y[!nas][ord]
    len <- length(dis[!nas])
    if(drop) {
        y <- y[-1]
        len <- len - 1
    }
    K <- seq_len(kmax)
    cumsum(y[K]) / K
}

##' @title Return the non-zero minimum of a vector of distances
##'
##' @param x the vector of distances for which the non-zero minimum is
##' required
##' @param drop logical; should the trivial (zero) distance of site
##' with itself be dropped?
##'
##' @return The minimum, non-zero distance, a vector of length 1.
##'
##' @author Gavin L. Simpson
minDij <- function(x, drop = TRUE) {
    ord <- order(x)
    if(drop)
        x[ord][2] # we don't want the first zero distance
    else
        x[ord][1]
}

##' @title The maximum bias statistic of transfer function residuals
##'
##' @param error numeric vector of model residuals
##' @param y numeric vector of observed environmental data
##' @param n numeric; number of sections to cut environmental gradient into
##'
##' @return A numeric vector of length \code{n} containing the value of
##' the largest residual in each of the \code{n} setions of the gradient.
##'
##' @author Gavin L. Simpson
maxBias <- function(error, y, n = 10) {
    groups <- cut.default(y, breaks = n, labels = 1:n)
    bias <- tapply(error, groups, mean)
    bias[which.max(abs(bias))]
}

##' @title Simple capitalisation function from ?toupper
##'
##' @param x string to be capitalised
##'
##' @return The capitalise string
##'
##' @author Gavin L. Simpson
.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep = "",
        collapse = " ")
}

###########################################################################
##                                                                       ##
## wmean - simple, quick version of weighted.mean                        ##
##                                                                       ##
## Created       : 16-Feb-2007                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 16-Feb-2007                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x - string to be capitalised                                          ##
##                                                                       ##
###########################################################################
wmean <- function(spp, env) {
  sum(env * spp)/sum(spp)
}

## w.avg - fast weighted mean function with no checks
`w.avg` <- function(x, env) {
    opt <- ColSums(x * env) / ColSums(x)
    names(opt) <- colnames(x)
    opt
}

## fast rowSums and colSums functions without the checking
`RowSums` <- function(x, na.rm = FALSE) {
    dn <- dim(x)
    p <- dn[2]
    dn <- dn[1]
    .rowSums(x, dn, p, na.rm)
}

`ColSums` <- function(x, na.rm = FALSE) {
    dn <- dim(x)
    n <- dn[1]
    dn <- dn[2]
    .colSums(x, n, dn, na.rm)
}

## w.tol --- computes weighted standard deviations AKA tolerances
w.tol <- function(x, env, opt, useN2 = TRUE) {
    ## x   = species abundances
    ## env = vector of response var
    ## opt = weighted average optima
    nr <- NROW(x)
    nc <- NCOL(x)
    tol <- .C("WTOL", x = as.double(env), w = as.double(x),
              opt = as.double(opt),
              nr = as.integer(nr), nc = as.integer(nc),
              stat = double(nc), NAOK = FALSE,
              PACKAGE="analogue")$stat
    if(useN2)
        tol <- tol / sqrt(1 - (1 / sppN2(x)))
    names(tol) <- colnames(x)
    tol
}

`sppN2` <- function(x) {
    ## quickly compute Hill's N2 for species
    ## x = species abundances
    ## ONLY call within an existing function
    tot <- ColSums(x)
    x <- sweep(x, 2, tot, "/")
    x <- x^2
    H <- ColSums(x, na.rm = TRUE)
    1/H
}

