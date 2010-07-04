## Fit a Principal Curve to matrix X
## Wrapper to principal.curve() in package princurve
## We use the original code plus our wrappers as pcurve()
## in package pcurve is too complex for our needs

## prcurve (named after prcomp): fits a principal curve to matrix X
prcurve <- function(X, method = c("ca","pca","random"),
                    smoother = smoothSpline,
                    complexity, vary = FALSE, maxComp,
                    finalCV = FALSE,
                    axis = 1, rank = FALSE, stretch = 2,
                    maxit = 10, trace = FALSE, thresh = 0.001,
                    plotit = FALSE, ...) {
    ## X should be a matrix, attempt to coerce
    if(!isTRUE(all.equal(class(X), "matrix")))
        X <- data.matrix(X)
    ## set/select default method for starting configuration
    if(missing(method))
        method <- "ca"
    else
        method <- match.arg(method)
    ## data stats
    n <- NROW(X) ## number of observations
    m <- NCOL(X) ## number of variables
    ## starting configuration
    config <- startConfig <- initCurve(X, method = method,
                                       rank = rank,
                                       axis = axis)
    ## Need to sort out auto DF choices after pcurve::pcurve
    ## Vary degrees of freedom per variable?
    if(missing(complexity)) {
        complexity <- numeric(length = m)
        for(j in seq_along(complexity)) {
            complexity[j] <- smoother(config$lambda, X[, j],
                                      choose = TRUE, ...)$complexity
        }
        if(!vary) {
            complexity <- rep(median(complexity), m)
        }
    } else {
        if((len <- length(complexity)) == 1) {
            complexity <- rep(complexity, m)
        } else if(len != m) {
            stop("Ambiguous 'complexity'; should be length 1 or NCOL(X)")
        }
    }
    if(missing(maxComp))
        maxComp <- 5 * log10(n)
    ## fix-upreset complexity > maxComp to maxComp
    complexity[complexity > maxComp] <- maxComp
    ##
    iter <- 0L
    if(trace)
        writeLines(sprintf("Initial curve: d.sq: %.4f", config$dist))
    ##dist.raw <- sum(diag(var(X))) * (NROW(X) - 1)
    dist.old <- sum(diag(var(X)))
    s <- matrix(NA, nrow = n, ncol = m)
    converged <- (abs((dist.old - config$dist)/dist.old) <=
                  thresh)
    ## Start iterations ----------------------------------------------
    while (!converged && iter < maxit) {
        iter <- iter + 1L
        for(j in seq_len(m)) {
            s[, j] <- fitted(smoother(config$lambda, X[, j],
                                      complexity = complexity[j],
                                      choose = FALSE, ...))
        }
        dist.old <- config$dist
        config <- get.lam(X, s = s, stretch = stretch)
        class(config) <- "prcurve"
        ## Converged?
        converged <- (abs((dist.old - config$dist)/dist.old) <=
                      thresh)
        if(plotit)
            plot(config, X, sub = paste("Iteration:", iter))
        if (trace)
            writeLines(sprintf(paste("Iteration %",
                                     max(3, nchar(maxit)),
                                     "i: d.sq: %.4f", sep = ""),
                               iter, config$dist))
    }
    ## End iterations ------------------------------------------------
    ## if we want a final CV spline fit?
    if(finalCV) {
        iter <- iter + 1L
        for(j in seq_len(n)) {
            sFit <- smoother(config$lambda, X[, j],
                             cv = TRUE, choose = TRUE, ...)
            s[, j] <- if(sFit$complexity > maxComp) {
                ## too complex, turn of CV and refit with max df allowed
                fitted(smoother(config$lambda, X[, j], cv = FALSE,
                                choose = FALSE,
                                complexity = maxComp,
                                ...))
            } else {
                fitted(sFit)
            }
        }
        config <- get.lam(X, s = config$s, stretch = stretch)
        class(config) <- "prcurve"
        if(plotit) {
            ## plot the iteration
            plot(config, X)
        }
        if (trace)
            writeLines(sprintf(paste("Iteration %",
                                     max(3, nchar(maxit)),
                                     "s: d.sq: %.4f", sep = ""),
                               "CV", config$dist))
    }
    names(config$tag) <- names(config$lambda) <-
        rownames(config$s) <- rownames(X)
    colnames(config$s) <- names(complexity) <- colnames(X)
    config$converged <- converged
    config$iter <- iter
    config$totalDist <- startConfig$dist
    config$complexity <- complexity
    config$call <- match.call()
    class(config) <- c("prcurve")
    return(config)
}
