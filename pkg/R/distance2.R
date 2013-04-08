## `distance2` <- function(x, ...)
##     UseMethod("distance2")

## `distance2.default` <- function(x, y,
##                                 method = c("euclidean",
##                                 "SQeuclidean", "chord",
##                                 "SQchord", "bray", "chi.square",
##                                 "SQchi.square", "information",
##                                 "chi.distance", "manhattan",
##                                 "kendall", "gower", "alt.gower",
##                                 "mixed"),
##                                 weights = NULL, R = NULL,
##                                 type = list(),
##                                 ordinal = c("gower","rank","metric"), ...) {
##     pColl <- function(n) paste(n, collapse = ", ")
##     ## Euclid?an could be spelled variously
##     if(!is.na(pmatch(method, "euclidian")))
## 	method <- "euclidean"
##     if (missing(method))
##         method <- "euclidean"
##     METHODS <- c("euclidean", "SQeuclidean", "chord", "SQchord",
##                  "bray", "chi.square", "SQchi.square",
##                  "information","chi.distance", "manhattan",
##                  "kendall", "gower", "alt.gower", "mixed")
##     method <- match.arg(method)
##     DCOEF <- pmatch(method, METHODS)
##     ordinal <- match.arg(ordinal)
##     ORDTYPES <- c("gower","rank","metric")
##     if(missing(y)) { ## only a single matrix
##         ## TODO
##     } else { ## two matrices
##         ## check x and y have same columns
##         if(!isTRUE(all.equal(names(x), names(y))))
##             stop("'x' and 'y' appear to have different variables.")
##         if(!isTRUE(all.equal((n.vars <- ncol(x)), ncol(y))))
##             stop("'x' and 'y' have different numbers of columns.")
##         ## variables
##         nrx <- nrow(x)
##         nry <- nrow(y)
##         d <- numeric(length = nrx * nry)
##         ## object names (row names)
##         x.names <- rownames(x)
##         y.names <- rownames(y)
##         ## some preprocessing steps required for some coefs
##         ## so dealt with separately
##         if(method %in% c("chi.distance", "gower", "alt.gower",
##                          "mixed", "kendall")) {
##             if(method == "chi.distance") {
##                 x <- data.matrix(x)
##                 y <- data.matrix(y)
##                 csum <- colSums(rbind(x, y))
##                 y <- y / rowSums(y)
##                 x <- x / rowSums(x)
##                 d <- .C("xy_chisq_dist", x = as.double(x), y = as.double(y),
##                         nr1 = as.integer(nrx), nr2 = as.integer(nry),
##                         nc = as.integer(n.vars), d = as.double(d),
##                         csum = as.double(csum), NAOK = as.integer(FALSE),
##                         PACKAGE = "analogue")$d
##             }
##             if(method %in% c("gower", "alt.gower", "mixed")) {
##                 if(method == "mixed") {
##                     if(is.null(weights))
##                         weights <- rep(1, n.vars)
##                     else {
##                         if(length(weights) != n.vars)
##                             stop("'weights' must be of length 'ncol(x)'")
##                     }
##                     ## process vtypes
##                     if(length(type)) {
##                         ## if 'type's supplied, validate
##                     }
##                     ## TODO
##                     if(is.data.frame(x)) {
##                         type2x <- sapply(x, data.class, USE.NAMES = FALSE)
##                         ##x <- data.matrix(x)
##                     } else {
##                         type2x <- rep("numeric", n.vars)
##                         names(type2x) <- colnames(x)
##                     }
##                     if(is.data.frame(y)) {
##                         type2y <- sapply(y, data.class, USE.NAMES = FALSE)
##                         ##y <- data.matrix(y)
##                     } else {
##                         type2y <- rep("numeric", n.vars)
##                         names(type2y) <- colnames(y)
##                     }
##                     ## x and y should have same column types
##                     if(!isTRUE(all.equal(type2x, type2y)))
##                         stop("Variable types in 'x' and 'y' differ.
## Did you forget  to 'join' 'x' and 'y' before calling 'distance'?")
##                     type2x[tI <- type2x %in% c("numeric", "integer")] <- "Q"
##                     ## save which are ordinal for rank conversion below
##                     type2x[(ordinal <- type2x == "ordered")] <- "O"
##                     type2x[type2x == "factor"] <- "N"
##                     type2x[type2x == "logical"] <- "A"
##                     typeCodes <- c("A", "S", "N", "O", "Q", "I", "T")
##                     type3 <- match(type2x, typeCodes)
##                     if (any(ina <- is.na(type3)))
##                         stop("invalid type ", type2x[ina], " for column numbers ",
##                              pColl(which(ina)))

##                     ## Convert to matrices from now on
##                     ## also takes care of ordinal == metric as all factors
##                     ## are converted to internal numeric codes
##                     x <- data.matrix(x)
##                     y <- data.matrix(y)

##                     ## Convert ordinal variables to ranks or numerics
##                     ## implemented as per Podani 1999. Only do ranks here as
##                     ## conversion to matrices above handled the standard case
##                     x[, ordinal] <- apply(x[, ordinal], 2, rank, na.last = "keep")
##                     y[, ordinal] <- apply(y[, ordinal], 2, rank, na.last = "keep")

##                     ## Compute range Rj
##                     XY <- rbind(x, y)
##                     if(is.null(R)) {
##                         maxi <- apply(XY, 2, max, na.rm = TRUE)
##                         mini <- apply(XY, 2, min, na.rm = TRUE)
##                         R <- maxi - mini
##                     } else {
##                         if(length(R) != n.vars)
##                             stop("'R' must be of length 'ncol(x)'")
##                     }

##                     ## For Ordinal we need TiMin and TiMax
##                     ## compute over all variables so they have same length as
##                     ## everything else
##                     doT <- function(X, which) {
##                         val <- if(which == "min") {
##                             min(X, na.rm = TRUE)
##                         } else {
##                             max(X, na.rm = TRUE)
##                         }
##                         nas <- is.na(X)
##                         length(which(X[!nas] == val))
##                     }
##                     tmin <- apply(XY, 2, doT, which = "min")
##                     tmax <- apply(XY, 2, doT, which = "max")

##                     ## How do we want to handle ordinals - convert to interger code
##                     ## for use in C
##                     podani <- match(ordinal, ORDTYPES)

##                     ## call the C code
##                     d <- .C("xy_mixed", x = as.double(x), y = as.double(y),
##                             nr1 = as.integer(nrx), nr2 = as.integer(nry),
##                             nc = as.integer(n.vars), d = as.double(d),
##                             vtype = as.integer(type3),
##                             weights = as.double(weights), R = as.double(R),
##                             tmin = as.integer(tmin), tmax = as.integer(tmax),
##                             podani = as.integer(podani),
##                             NAOK = as.integer(TRUE),
##                             PACKAGE = "analogue")$d
##                 } else {
##                     if(is.null(R)) {
##                         XY <- rbind(x, y)
##                         maxi <- apply(XY, 2, max, na.rm = TRUE)
##                         mini <- apply(XY, 2, min, na.rm = TRUE)
##                         R <- maxi - mini
##                     } else {
##                         if(length(R) != n.vars)
##                             stop("'R' must be of length 'ncol(x)'")
##                     }
##                     x <- data.matrix(x)
##                     y <- data.matrix(y)

##                     ## pre-process for gower and alt gower
##                     ## but these handled by xy_distance below
##                     x <- sweep(x, 2, R, "/")
##                     y <- sweep(y, 2, R, "/")
##                     d <- .C("xy_distance", x = as.double(x), y = as.double(y),
##                             nr1 = as.integer(nrx), nr2 = as.integer(nry),
##                             nc = as.integer(n.vars), d = as.double(d),
##                             method = as.integer(DCOEF), NAOK = as.integer(FALSE),
##                             PACKAGE = "analogue")$d
##                 }
##             }
##             if(method == "kendall") {
##                 x <- data.matrix(x)
##                 y <- data.matrix(y)
##                 XY <- rbind(x, y)
##                 maxi <- apply(XY, 2, max)
##                 d <- .C("xy_kendall", x = as.double(x), y = as.double(y),
##                         nr1 = as.integer(nrx), nr2 = as.integer(nry),
##                         nc = as.integer(n.vars), d = as.double(d),
##                         maxi = as.double(maxi), NAOK = as.integer(FALSE),
##                         PACKAGE = "analogue")$d
##             }
##         } else {
##             ## must be one of the DC's handled by xy_distance
##             x <- data.matrix(x)
##             y <- data.matrix(y)
##             d <- .C("xy_distance", x = as.double(x), y = as.double(y),
##                     nr1 = as.integer(nrx), nr2 = as.integer(nry),
##                     nc = as.integer(n.vars), d = as.double(d),
##                     method = as.integer(DCOEF), NAOK = as.integer(FALSE),
##                     PACKAGE = "analogue")$d
##         }
##         ## convert d to a matrix
##         d <- matrix(d, ncol = n.vars, byrow = TRUE)
##         colnames(d) <- y.names
##         rownames(d) <- x.names
##         attr(d, "method") <- method
##         class(d) <- c("distance","matrix")
##     }
##     return(d)
## }

## set.seed(1)
## bar <- matrix(sample(3, 9, replace = TRUE), ncol = 3)
## foo <- matrix(sample(3, 9, replace = TRUE), ncol = 3)
## foobar <- rbind(bar, foo)
## out <- matrix(ncol = ncol(bar), nrow = nrow(foobar))
## res <- numeric(length = nrow(foobar))
## for(i in seq_len(nrow(foobar))) {
##     for(j in seq_len(ncol(foobar))) {
##         for(k in seq_along(foobar[,j])) {
##             res[k] <- foobar[k,j] == foobar[i,j]
##         }
##         out[i, j] <- sum(res)
##     }
## }

## set.seed(1)
## bar <- matrix(sample(3, 9, replace = TRUE), ncol = 3)
## foo <- matrix(sample(3, 9, replace = TRUE), ncol = 3)
## outbar <- matrix(0, ncol = ncol(bar), nrow = nrow(bar))
## outfoo <- matrix(0, ncol = ncol(foo), nrow = nrow(foo))
## resbar <- numeric(length = nrow(bar))# + nrow(foo))
## resfoo <- numeric(length = nrow(bar))# + nrow(foo))

## for(i in seq_len(ncol(bar))) {
##     for(j in seq_len(nrow(bar))) {
##         for(k in seq_len(nrow(bar))) {
##             resbar[k] <- bar[k, i] == bar[j, i]
##         }
##         outbar[j, i] <- sum(resbar)
##     }
##     for(j in seq_len(nrow(foo))) {
##         for(k in seq_len(nrow(foo))) {
##             resfoo[k] <- foo[k, i] == bar[j, i]
##         }
##         outfoo[j, i] <- sum(resfoo)
##     }
## }

`distance3` <- function(x, ...)
    UseMethod("distance3")

`distance3.default` <-
    function(x, y, method = "euclidean",
             weights = NULL, R = NULL,
             ...)
{
    METHODS <- c("euclidean", "SQeuclidean", "chord", "SQchord",
                 "bray", "chi.square", "SQchi.square", "information",
                 "chi.distance", "manhattan", "kendall", "gower", "alt.gower",
                 "mixed")
    pColl <- function(n) paste(n, collapse = ", ")
    ## Euclid?an could be spelled variously
    if(!is.na(pmatch(method, "euclidian")))
	method <- "euclidean"
    METHODS <- c("euclidean", "SQeuclidean", "chord", "SQchord",
                 "bray", "chi.square", "SQchi.square",
                 "information","chi.distance", "manhattan",
                 "kendall", "gower", "alt.gower", "mixed")
    DCOEF <- pmatch(method, METHODS)
    if(missing(y)) { ## only a single matrix
        ## TODO
    } else { ## two matrices
        ## check x and y have same columns
        if(!isTRUE(all.equal(names(x), names(y))))
            stop("'x' and 'y' appear to have different variables.")
        if(!isTRUE(all.equal((n.vars <- ncol(x)), ncol(y))))
            stop("'x' and 'y' have different numbers of columns.")
        ## variables
        nrx <- nrow(x)
        nry <- nrow(y)
        d <- numeric(length = nrx * nry)
        ## object names (row names)
        x.names <- rownames(x)
        y.names <- rownames(y)
        ## some preprocessing steps required for some coefs
        ## so dealt with separately
        if(method %in% c("chi.distance", "gower", "alt.gower",
                         "mixed", "kendall")) {
            if(method == "chi.distance") {
                x <- data.matrix(x)
                y <- data.matrix(y)
                csum <- colSums(rbind(x, y))
                y <- y / rowSums(y)
                x <- x / rowSums(x)
                d <- .C("xy_chisq_dist", x = as.double(x), y = as.double(y),
                        nr1 = as.integer(nrx), nr2 = as.integer(nry),
                        nc = as.integer(n.vars), d = as.double(d),
                        csum = as.double(csum), NAOK = as.integer(FALSE),
                        PACKAGE = "analogue")$d
            }
            if(method %in% c("gower", "alt.gower", "mixed")) {
                if(method == "mixed") {
                    if(is.null(weights))
                        weights <- rep(1, n.vars)
                    else {
                        if(length(weights) != n.vars)
                            stop("'weights' must be of length 'ncol(x)'")
                    }
                    ## process vtypes
                    if(is.data.frame(x)) {
                        type2x <- sapply(x, data.class, USE.NAMES = FALSE)
                        ##x <- data.matrix(x)
                    } else {
                        type2x <- rep("numeric", n.vars)
                        names(type2x) <- colnames(x)
                    }
                    if(is.data.frame(y)) {
                        type2y <- sapply(y, data.class, USE.NAMES = FALSE)
                        ##y <- data.matrix(y)
                    } else {
                        type2y <- rep("numeric", n.vars)
                        names(type2y) <- colnames(y)
                    }
                    ## x and y should have same column types
                    if(!isTRUE(all.equal(type2x, type2y)))
                        stop("Variable types in 'x' and 'y' differ.
Did you forget  to 'join' 'x' and 'y' before calling 'distance'?")

                    ## Record the variable types
                    type2x[tI <- type2x %in% c("numeric", "integer")] <- "Q"
                    ## save which are ordinal for rank conversion below - TODO
                    type2x[(ordinal <- type2x == "ordered")] <- "O"
                    type2x[type2x == "factor"] <- "N"
                    type2x[type2x == "logical"] <- "A"
                    typeCodes <- c("A", "S", "N", "O", "Q", "I", "T")
                    type3 <- match(type2x, typeCodes)
                    if (any(ina <- is.na(type3)))
                        stop("invalid type ", type2x[ina], " for column numbers ",
                             pColl(which(ina)))

                    ## Convert to matrices from now on
                    ## also takes care of ordinal == metric as all factors
                    ## are converted to internal numeric codes
                    x <- data.matrix(x)
                    y <- data.matrix(y)

                    ## Compute range Rj
                    XY <- rbind(x, y)
                    if(is.null(R)) {
                        maxi <- apply(XY, 2, max, na.rm = TRUE)
                        mini <- apply(XY, 2, min, na.rm = TRUE)
                        R <- maxi - mini
                    } else {
                        if(length(R) != n.vars)
                            stop("'R' must be of length 'ncol(x)'")
                    }

                    ## call the C code
                    d <- .C("xy_mixed", x = as.double(x), y = as.double(y),
                            nr1 = as.integer(nrx), nr2 = as.integer(nry),
                            nc = as.integer(n.vars), d = as.double(d),
                            vtype = as.integer(type3),
                            weights = as.double(weights), R = as.double(R),
                            NAOK = as.integer(TRUE),
                            PACKAGE = "analogue")$d
                } else {
                    if(is.null(R)) {
                        XY <- rbind(x, y)
                        maxi <- apply(XY, 2, max, na.rm = TRUE)
                        mini <- apply(XY, 2, min, na.rm = TRUE)
                        R <- maxi - mini
                    } else {
                        if(length(R) != n.vars)
                            stop("'R' must be of length 'ncol(x)'")
                    }
                    x <- data.matrix(x)
                    y <- data.matrix(y)

                    ## pre-process for gower and alt gower
                    ## but these handled by xy_distance below
                    x <- sweep(x, 2, R, "/")
                    y <- sweep(y, 2, R, "/")
                    d <- .C("xy_distance", x = as.double(x), y = as.double(y),
                            nr1 = as.integer(nrx), nr2 = as.integer(nry),
                            nc = as.integer(n.vars), d = as.double(d),
                            method = as.integer(DCOEF), NAOK = as.integer(FALSE),
                            PACKAGE = "analogue")$d
                }
            }
            if(method == "kendall") {
                x <- data.matrix(x)
                y <- data.matrix(y)
                XY <- rbind(x, y)
                maxi <- apply(XY, 2, max)
                d <- .C("xy_kendall", x = as.double(x), y = as.double(y),
                        nr1 = as.integer(nrx), nr2 = as.integer(nry),
                        nc = as.integer(n.vars), d = as.double(d),
                        maxi = as.double(maxi), NAOK = as.integer(FALSE),
                        PACKAGE = "analogue")$d
            }
        } else {
            ## must be one of the DC's handled by xy_distance
            x <- data.matrix(x)
            y <- data.matrix(y)
            d <- .C("xy_distance", x = as.double(x), y = as.double(y),
                    nr1 = as.integer(nrx), nr2 = as.integer(nry),
                    nc = as.integer(n.vars), d = as.double(d),
                    method = as.integer(DCOEF), NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
        ## convert d to a matrix
        d <- matrix(d, ncol = nry, byrow = TRUE)
        colnames(d) <- y.names
        rownames(d) <- x.names
        attr(d, "method") <- method
        attr(d, "type") <- "asymmetric"
        class(d) <- c("distance","matrix")
    }
    return(d)
}
