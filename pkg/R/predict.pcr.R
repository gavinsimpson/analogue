`predict.pcr` <- function(object, newdata, ncomp = seq_len(object$ncomp),
                          CV = c("none", "LOO", "bootstrap", "nfold"),
                          verbose = FALSE, n.boot = 100, nfold = 5,
                          ...) {
    if(missing(newdata))
        return(fitted(object))
    ## store names of new samples
    newSamp <- rownames(newdata)
    newdata <- as.matrix(newdata)
    if (missing(CV))
        CV <- "none"
    CV <- match.arg(CV)
    Np <- NROW(newdata)
    B <- coef(object)
    if(identical(CV, "none")) {
        want <- (spp.names <- colnames(object$data$x)) %in% colnames(newdata)
        want <- spp.names[want]
        newdata <- newdata[, want, drop = FALSE]
        ## apply transformation to newdata
        newdata <- object$tranFun(newdata)
        ## do predictions
        ## matrix of predictions
        pred <- matrix(ncol = length(ncomp), nrow = Np)
        for(j in ncomp) {
            B0 <- object$yMean - object$xMeans %*% B[, j]
            pred[, j] <- newdata %*% B[, j] + rep(B0, Np)
        }
    } else {
        stop("Other methods of crossvalidation not yet implemented")
    }
    rownames(pred) <- newSamp
    colnames(pred) <- paste0("PC", ncomp)
    pred
}

