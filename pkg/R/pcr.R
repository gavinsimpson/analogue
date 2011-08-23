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
    fun.supplied <- FALSE
    if(!missing(tranFun)) {
        FUN <- match.fun(tranFun)
        x <- tranFun(x)
        tranParms <- attr(x, "parms")
        fun.supplied <- TRUE
    } else {
        FUN <- tranParms <- NA
        tranFun <- "none"
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
                              RMSE = sqrt(colMeans(residuals^2)))

    ## apply some names to prettify to output
    rownames(performance) <- colnames(B) <- colnames(fitted.values) <-
        colnames(residuals) <- names(varExpl) <- colnames(TT) <-
            colnames(P) <- rownames(tQ) <- paste("PC", S, sep = "")
    rownames(B) <- rownames(P) <- colnames(x)
    rownames(fitted.values) <- rownames(residuals) <- rownames(TT) <- rownames(x)

    ## get and fix up the call
    .call <- match.call()
    .call[[1]] <- as.name("pcr")
    if(fun.supplied) {
        ## fix-up the name of the transformation function used,
        ## needed when formula method called...
        .call[[which(names(.call) == "tranFun")]] <- as.name(deparse(substitute(tranFun)))
    }

    ## return object
    Obj <- list(fitted.values = fitted.values,
                coefficients = B,
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
    .call <- match.call()
    .call[[1]] <- as.name("pcr")
    ## fix-up the name of the transformation function used,
    ## needed when formula method called...
    ##.call[[which(names(.call) == "tranFun")]] <- as.name(deparse(substitute(tranFun)))
    Obj$call <- .call
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
    writeLines("RMSE (Apparent):")
    perf <- x$performance[, "RMSE"]
    names(perf) <- rownames(x$performance)
    print(perf)
}

`screeplot.pcr` <- function(x, restrict = NULL,
                            display = c("RMSE","avgBias","maxBias","R2"),
                            xlab = NULL, ylab = NULL, main = NULL, sub = NULL,
                            ...) {
    display <- match.arg(display)
    captions <- c("RMSE", "Average bias", "Maximum bias", "R squared")
    names(captions) <- c("RMSE", "avgBias", "maxBias", "R2")
    if (is.null(xlab))
        xlab <- "No. of components"
    if (is.null(ylab))
        ylab <- captions[display]
    if (is.null(main))
        main <- deparse(substitute(x))
    if (is.null(sub)) {
        cal <- x$call
        ##if (!is.na(m.f <- match("formula", names(cal)))) {
        ##    cal <- cal[c(1, m.f)]
        ##    names(cal)[2] <- ""
        ##}
        cc <- deparse(cal, 90)
        nc <- nchar(cc[1])
        abbr <- length(cc) > 1 || nc > 90
        sub <- if (abbr)
            paste(substr(cc[1], 1, min(90, nc)), "...")
        else cc[1]
    }
    dat <- performance(x)
    if(!is.null(restrict)) {
        comps <- min(restrict, x$ncomp)
    } else {
        comps <- x$ncomp
    }
    Scomps <- seq_len(comps)
    plot(Scomps, dat[Scomps, display], type = "n", ylab = ylab, xlab = xlab,
         main = main, sub = sub, ...)
    if(comps > 20) {
        lines(Scomps, dat[Scomps, display], type = "b", ...)
    } else {
        lines(Scomps, dat[Scomps, display], type = "b", ..., pch = NA)
        text(Scomps, dat[Scomps, display], labels = as.character(Scomps), cex = 0.8,
             ...)
    }
    invisible()
}

`performance.pcr` <- function(object, ...) {
    retval <- object$performance
    class(retval) <- c("performance","data.frame")
    retval
}

`coef.pcr` <- function(object, comps = NULL, ...) {
    coefs <- object$coefficients
    nc <- NCOL(coefs)
    if(is.null(comps))
        comps <- seq_len(nc)
    else {
        if(!is.numeric(comps))
            stop("Non-numeric selection of components requested.")
        if(min(comps) < 1 || max(comps) > nc)
            stop("Requested components outside range of actual components.")
    }
    coefs <- coefs[, comps]
    coefs
}

`fitted.pcr` <- function(object, comps = NULL, ...) {
    fits <- object$fitted.values
    nc <- NCOL(fits)
    if(is.null(comps))
        comps <- seq_len(nc)
    else {
        if(!is.numeric(comps))
            stop("Non-numeric selection of components requested.")
        if(min(comps) < 1 || max(comps) > nc)
            stop("Requested components outside range of actual components.")
    }
    fits <- fits[, comps]
    fits
}

`residuals.pcr` <- function(object, comps = NULL, ...) {
    resi <- object$residuals
    nc <- NCOL(resi)
    if(is.null(comps))
        comps <- seq_len(nc)
    else {
        if(!is.numeric(comps))
            stop("Non-numeric selection of components requested.")
        if(min(comps) < 1 || max(comps) > nc)
            stop("Requested components outside range of actual components.")
    }
    resi <- resi[, comps]
    resi
}

`eigenvals.pcr` <- function(x, ...) {
    x$varExpl
}
