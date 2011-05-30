## Principal Components Regression Models

## generic
`pcr` <- function(x, ...) {
    UseMethod("pcr")
}

## default
`pcr.default` <- function(x, y, ncomp, tranFun, ...) {
    ## convert to matrices for speed
    x <- data.matrix(x)

    ## dimensions
    Nx <- NROW(x)
    Ny <- length(y)
    Mx <- NCOL(x) ## number of predictor vars

    ## Apply a transformation if specified
    if(!missing(tranFun)) {
        FUN <- match.fun(tranFun)
        x <- tranFun(x)
        tranParms <- attr(x, "parms")
    } else {
        FUN <- NA
    }

    ## centre x and y
    xMeans <- colMeans(x)
    yMean <- mean(y)
    x <- sweep(x, 2, xMeans, "-")
    y <- y - yMean

    ## How many components?
    ncomp <- if(missing(ncomp)) {
        min(Nx - 1, Mx)
    } else {
        if(ncomp < 1 || ncomp > (newcomp <- min(Nx - 1, Mx)))
            warning("Invalid 'ncomp'. Resetting to max possible.")
        newcomp
    }
    S <- seq_len(ncomp)

    ## Storage objects:
    ## B model coefficients
    ## fitted.values
    B <- matrix(0, nrow = Mx, ncol = ncomp)
    fitted.values <- matrix(0, nrow = Nx, ncol = ncomp)

    ## SVD
    SVD <- La.svd(x)
    D <- SVD$d[S]
    TT <- SVD$u[, S, drop = FALSE] %*% diag(D, nrow = ncomp)
    P <- t(SVD$vt[S, , drop = FALSE])
    tQ <- crossprod(TT, y) / (varExpl <- D^2)

    ## compute coefficients
    for(b in S) {
        bS <- seq_len(b)
        B[, b] <- P[, bS, drop = FALSE] %*% tQ[bS, ]
        fitted.values[, b] <- TT[, bS, drop = FALSE] %*% tQ[bS, ]
    }

    ## other model output
    residuals <- y - fitted.values
    fitted.values <- sweep(fitted.values, 2, yMean, "+")

    ## model performance
    Y <- y + yMean
    performance <- data.frame(R2 = drop(cor(fitted.values, Y)),
                              avgBias = colMeans(residuals),
                              maxBias = apply(residuals, 2, maxBias, Y),
                              RMSEP = sqrt(colMeans(residuals^2)))
    rownames(performance) <- paste("PC", S, sep = "")

    ## get and fix up the call
    .call <- match.call()
    .call[[1]] <- as.name("pcr")
    ## fix-up the name of the transformation function used,
    ## needed when formula method called...
    .call[[which(names(.call) == "tranFun")]] <- as.name(deparse(substitute(tranFun)))

    ## return object
    Obj <- list(fitted.values = fitted.values,
                coefficients = coefficients,
                residuals = residuals,
                scores = TT,
                loadings = P,
                Yloadings = t(tQ),
                xMeans = xMeans,
                yMean = yMean,
                varExpl = varExpl,
                totvar = sum(x * x),
                call = .call,
                tranFun = FUN,
                tranParms = tranParms,
                performance = performance,
                ncomp = ncomp)
    class(Obj) <- "pcr"
    Obj
}

`pcr.formula` <- function(formula, data, subset, na.action,
                          ..., model = FALSE) {
    ## the function call
    .call <- match.call()
    ## need to reset due to method dispatch
    .call[[1]] <- as.name("pcr")

    ## keep only the arguments which should go into the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    attr(attr(mf, "terms"), "intercept") <- 0

    ## terms objects
    mt <- attr(mf, "terms")
    ## model matrices
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)

    ## fit & the PCR
    Obj <- pcr.default(x = x, y = y, ...)
    Obj$na.action <- attr(mf, "na.action")
    Obj$terms <- mt
    if(model)
        Obj$model <- mf
    Obj
}

`Hellinger` <- function(x, apply = FALSE) {
    tran(x, method = "hellinger")
    ## ignore 'apply' as no meta-parameters
}

`ChiSquare` <- function(x, apply = FALSE) {
    if(apply) {
        ## apply pre-computed meta-parameters to transform
        ## test samples to match training samples
        parms <- attr(x, "parms")
        x <- with(parms, sqrt(gsum) * x/outer(rsum, sqrt(csum)))
    } else {
        ## perform transformation and preserve the meta-parameters
        if (any(x < 0, na.rm = TRUE)) {
            k <- min(x, na.rm = TRUE)
            warning("'x'contains negative entries: result may be nonsense\n")
        } else {
            k <- .Machine$double.eps
        }
        gsum <- sum(x, na.rm = TRUE)
        rsum <- pmax(k, rowSums(x, na.rm = TRUE))
        csum <- colSums(x, na.rm = TRUE)
        x <- sqrt(gsum) * x/outer(rsum, sqrt(csum))
        attr(x, "parms") <- list(gsum = gsum, rsum = rsum, csum = csum)
    }
    x
}

`print.pcr` <- function(x, digits = min(getOption("digits"), 4), ...) {
    cat("\n")
    writeLines(strwrap("Principal Component Regression Model", prefix = "\t"),
               sep = "\n\n")
    writeLines(strwrap("Call:"))
    print(x$call)
    cat("\n")
    writeLines(strwrap(paste("No. of Components:", x$ncomp)), sep = "\n\n")
    writeLines("RMSEP (Apparent):")
    perf <- x$performance[, "RMSEP"]
    names(perf) <- rownames(x$performance)
    print(perf)
}
