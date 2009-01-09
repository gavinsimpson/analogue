###########################################################################
##                                                                       ##
## distance - function to compute distances between samples              ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 13-Oct-2007                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
## x = training data, y = fossil data
distance <- function(x, ...) UseMethod("distance")

distance.join <- function(x, ...)
  {
    if(!inherits(x, "join"))
      stop("This method should only be used on objects of class 'join'")
    if(inherits(x, "data.frame")) {
      distance.default(x, ...)
    } else {
      if(length(x) != 2)
        warning("Object contains more than 2 data sets.\n  Only the first 2 data sets used")
      distance.default(x[[1]], x[[2]], ...)
    }
  }

distance.default <- function(x, y,
                             method = c("euclidean", "SQeuclidean", "chord",
                               "SQchord", "bray", "chi.square", "SQchi.square",
                               "information", "chi.distance", "manhattan",
                               "kendall", "gower", "alt.gower", "mixed"),
                             fast = TRUE,
                             weights = NULL, R = NULL, ...)
  {
    euclidean <- function(x, y)
      {
        sqrt(sum((x - y)^2))
      }
    SQeuclidean <- function(x, y)
      {
        sum((x - y)^2)
      }
    chord <- function(x, y)
      {
        x <- sqrt(x); y <- sqrt(y)
        euclidean(x, y)
      }
    SQchord <- function(x, y)
      {
        x <- sqrt(x); y <- sqrt(y)
        SQeuclidean(x, y)
      }
    bray <- function(x, y)
      {
        sum(abs(x - y)) / sum(x + y)
      }
    chi.square <- function(x, y)
      {
        inds <- !(x == 0 & y == 0)
        sqrt(sum(((x[inds] - y[inds])^2) / (x[inds] + y[inds])))
      }
    SQchi.square <- function(x, y)
      {
        inds <- !(x == 0 & y == 0)
        sum(((x[inds] - y[inds])^2) / (x[inds] + y[inds]))
      }
    information <- function(x, y)
      {
        XY <- x + y
        A <- x * log2((2 * x) / XY)
        B <- y * log2((2 * y) / XY)
        sum(A, na.rm = TRUE) + sum(B, na.rm = TRUE)
      }
    chi.distance <- function(x, y, colsum)
      {
        sqrt(sum(((x - y)^2) / (colsum / sum(colsum))))
      }
    manhattan <- function(x, y)
      {
        sum(abs(x - y))
      }
    kendall <- function(x, y, maxi)
      {
        sum(maxi - min(x, y))
      }
    gower <- function(x, y, maxi, mini)
      {
        sum(abs(x - y) / (maxi - mini), na.rm = TRUE)
      }
    alt.gower <- function(x, y, maxi, mini)
      {
        sqrt(2 * sum(abs(x - y) / (maxi - mini), na.rm = TRUE))
      }
    mixed <- function(x, y, facs, weights, R)
      {
        weights[is.na(x)] <- 0
        weights[is.na(y)] <- 0
        if(any(!facs)) {
          quant <- 1 - (abs(x[!facs] - y[!facs]) / R[!facs])
          quant <- quant * weights[!facs]
        }
        if(any(facs)) {
          factors <- ifelse(x[facs] == y[facs], 1, 0)
          factors <- factors * weights[facs]
        }
        if(any(!facs)) {
          if(any(facs))
            retval <- sum(factors, quant, na.rm = TRUE) / sum(weights)
          else
            retval <- sum(quant, na.rm = TRUE) / sum(weights)
        } else {
          retval <- sum(factors, na.rm = TRUE) / sum(weights)
        }
        return(1- retval)
      }
    Dist <- function(y, x, method, ...)#colsum = NULL)
      {
        dotargs <- list(...)
        switch(method,
               euclidean = apply(x, 1, euclidean, y),
               SQeuclidean = apply(x, 1, SQeuclidean, y),
               chord = apply(x, 1, chord, y),
               SQchord = apply(x, 1, SQchord, y),
               bray = apply(x, 1, bray, y),
               chi.square = apply(x, 1, chi.square, y),
               SQchi.square = apply(x, 1, SQchi.square, y),
               information = apply(x, 1, information, y),
               chi.distance = apply(x, 1, chi.distance, y,
                 colsum = dotargs$colsum),
               manhattan = apply(x, 1, manhattan, y),
               kendall = apply(x, 1, kendall, y, maxi = dotargs$maxi),
               gower = apply(x, 1, gower, y, maxi = dotargs$maxi,
                 mini = dotargs$mini),
               alt.gower = apply(x, 1, alt.gower, y, maxi = dotargs$maxi,
                 mini = dotargs$mini),
               mixed = apply(x, 1, mixed, y, facs = dotargs$facs,
                 weights = dotargs$weights, R = dotargs$R)
               )
      }
    if(missing(method))
      method <- "euclidean"
    method <- match.arg(method)
    y.miss <- FALSE
    if(missing(y)) {
      y.miss <- TRUE
      y <- x
    }
    if(method == "mixed") {
      ## sanity check: are same columns in x and y factors
      facs.x <- sapply(as.data.frame(x), is.factor)
      facs.y <- sapply(as.data.frame(y), is.factor)
      if(sum(facs.x - facs.y) > 0)
        stop("Different columns (species) are coded as factors in 'x' and 'y'")
      ## sanity check: levels of factors also need to be the same
      for(i in seq_along(facs.x)[facs.x]){
        if(!identical(levels(x[,i]), levels(y[,i])))
          stop("The levels of one or more factors in 'x' and 'y'\ndo not match.\nConsider using 'join(x, y)'. See '?join'")
      }
    }
    x.names <- rownames(x)
    x <- data.matrix(x)
    n.vars <- ncol(x)
    ## Do we want to remove NAs? Yes if gower, alt.gower and mixed,
    ## but fail for others
    NA.RM <- FALSE
    if(method %in% c("gower", "alt.gower", "mixed"))
      NA.RM <- TRUE
    #y.miss <- FALSE
    if(missing(y)) {
      #colsumx <- colSums(x, na.rm = NA.RM)
      #if(any(colsumx <= 0)) {
      #  x <- x[, colsumx > 0, drop = FALSE]
      #  warning("some species contain no data and were removed from data matrix x\n")
      #}
      y.miss <- TRUE
      y <- x
      y.names <- x.names
    } else {
      #if(method == "mixed") {
        ## sanity check: are same columns in x and y factors
        #facs.y <- sapply(as.data.frame(y), is.factor)
        #if(sum(facs.x - facs.y) > 0)
        #  stop("Different columns (species) are coded as factors in 'x' and 'y'")
        ## sanity check: levels of factors also need to be the same
        #for(i in seq_along(facs.x)[facs.x]){
        #  if(!identical(levels(x[,i]), levels(y[,i])))
        #    stop("The levels of one or more factors in 'x' and 'y' do not match.\nConsider using 'join(x, y)'. See '?join'")
        #}
      #}
      y.names <- rownames(y)
      y <- data.matrix(y)
    }
    if(method == "chi.distance")
      colsum <- colSums(join(as.data.frame(x),as.data.frame(y), split = FALSE))
    if(method == "mixed") {
      ## sort out the weights used, eg the Kroneker's Deltas
      ## weights must be NULL or numeric vector of length == ncol(x)
      if(is.null(weights))
        weights <- rep(1, n.vars)
      else {
        if(length(weights) != n.vars)
          stop("'weights' must be of length 'ncol(x)'")
      }
    }
    if(method == "kendall") {
        maxi <- apply(rbind(apply(x, 2, max), apply(y, 2, max)),
                      2, max)
    }
    if(method %in% c("gower", "alt.gower", "mixed")) {
        maxi <- mini <- numeric(length = n.vars)
        maxi <- apply(rbind(apply(x, 2, max),
                            apply(y, 2, max)),
                      2, max, na.rm = TRUE)
        mini <- apply(rbind(apply(x, 2, min),
                            apply(y, 2, min)),
                      2, min, na.rm = TRUE)
        if(is.null(R))
            R <- maxi - mini
        else {
            if(length(R) != n.vars)
                stop("'R' must be of length 'ncol(x)'")
        }
    }
    dimnames(x) <- dimnames(y) <- NULL
    if(method == "chi.distance") {
      y <- y / rowSums(y)
      x <- x / rowSums(x)
      res <- apply(y, 1, Dist, x, method, colsum = colsum)
    } else if (method == "kendall") {
      res <- apply(y, 1, Dist, x, method, maxi = maxi)
    } else if (method %in% c("gower", "alt.gower")) {
      res <- apply(y, 1, Dist, x, method, maxi = maxi, mini = mini)
    } else if (method == "mixed") {
      res <- apply(y, 1, Dist, x, method, facs = facs.x, weights = weights,
                   R = R)
    } else if (method %in% c("euclidean","SQeuclidean","chord","SQchord","bray") &&
               fast == TRUE && y.miss == TRUE) {
      if(method == "bray") {
        res <- as.matrix(vegdist(x, method = "bray"))
      } else {
        if(method %in% c("chord", "SQchord"))
          res <- as.matrix(dist(sqrt(x), method = "euclidean"))
        else
          res <- as.matrix(dist(x, method = "euclidean"))
        if(method %in% c("SQeuclidean", "SQchord"))
          res <- res^2
      }
    } else {
      res <- apply(y, 1, Dist, x, method)
    }
    if(is.null(dim(res))) {
      names(res) <- x.names
    } else {
      colnames(res) <- y.names
      rownames(res) <- x.names
    }
    attr(res, "method") <- method
    return(res)
  }
