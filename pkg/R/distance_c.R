distanceX <- function(x, y, method = "euclidean", weights = NULL, R = NULL,
                      as.dist = FALSE, ...) {
    ## Euclid?an could be spelled variously
    if(!is.na(pmatch(method, "euclidian")))
	method <- "euclidean"
    METHODS <- c("euclidean", "SQeuclidean", "chord", "SQchord",
                 "bray", "chi.square", "SQchi.square",
                 "information","chi.distance", "manhattan",
                 "kendall", "gower", "alt.gower", "mixed")
    DCOEF <- pmatch(method, METHODS)
    if(missing(y)) { ## only a single matrix
        ## variables
        nr <- nrow(x)
        nc <- ncol(x)
        ## object names (row names)
        x.names <- rownames(x)
        ## some preprocessing steps required for some coefs
        ## so dealt with separately
        if(method %in% c("chi.distance", "gower", "alt.gower",
                         "mixed", "kendall")) {
          if(method == "chi.distance") {
            x <- data.matrix(x)
            csum <- colSums(x)
            x <- x / rowSums(x)
            d <- .Call("Cchisqdistxx", x, csum, PACKAGE = "analogue")
          }
          if(method == "kendall") {
            x <- data.matrix(x)
            maxi <- apply(x, 2, max)
            d <- .Call("Ckendallxx", x, maxi, PACKAGE = "analogue")
          }
          if(method %in% c("gower", "alt.gower")) {
            if(is.null(R)) {
              x <- data.matrix(x)
              maxi <- apply(x, 2, max, na.rm = TRUE)
              mini <- apply(x, 2, min, na.rm = TRUE)
              R <- maxi - mini
            } else {
              if(length(R) != nc)
                stop("'R' must be of length 'ncol(x)'")
            }
            ## pre-process here for gower and alt.gower
            ## but note we call the main driver Cdistxx
            x <- sweep(x, 2, R, "/")
            d <- .Call("Cdistxx", x, DCOEF, PACKAGE = "analogue")
          }
        } else {
            ## must be one of the DC's handled by xy_distance
            x <- data.matrix(x)
            d <- .Call("Cdistxx", x, DCOEF, PACKAGE = "analogue")
        }
        attr(d, "Size") <- nr
        attr(d, "Labels") <- x.names
        attr(d, "Diag") <- FALSE
        attr(d, "Upper") <- FALSE
        attr(d, "method") <- method
        attr(d, "call") <- match.call()
        class(d) <- "dist"
        if(!as.dist) {
            d <- as.matrix(d)
            attr(d, "method") <- method
            attr(d, "type") <- "symmetric"
            class(d) <- c("distance","matrix")
        }
    } else { ## two matrices
        ## check x and y have same columns
        if(!isTRUE(all.equal(names(x), names(y))))
            stop("'x' and 'y' appear to have different variables.")
        if(!isTRUE(all.equal((n.vars <- ncol(x)), ncol(y))))
            stop("'x' and 'y' have different numbers of columns.")
        ## variables
        nrx <- nrow(x)
        nry <- nrow(y)
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
            ##d <- .Call("Cchisqdistxy", x, y, )
          }
        } else {
          ## must be one of the DC's handled by xy_distance
          x <- data.matrix(x)
          y <- data.matrix(y)
          d <- .Call("Cdistxy", x, y, DCOEF, PACKAGE = "analogue")
        }

        ## convert d to a matrix
        d <- matrix(d, ncol = nry, byrow = TRUE)
        colnames(d) <- y.names
        rownames(d) <- x.names
        attr(d, "method") <- method
        attr(d, "type") <- "asymmetric"
        class(d) <- c("distance","matrix")
    }
    d
}
