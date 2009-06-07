`predict.wa` <- function(object, newdata,
                         CV = c("none","LOO","bootstrap", "nfold"),
                         verbose = FALSE, n.boot = 100, nfold = 5,
                         ...) {
    if(missing(newdata))
        return(fitted(object))
    newdata <- as.matrix(newdata)
    if(missing(CV))
        CV <- "none"
    CV <- match.arg(CV)
    ## summaries
    spp.train <- object$n.spp
    spp.fossil <- ncol(newdata)
    n.train <- object$n.samp
    n.fossil <- nrow(newdata)
    n.in.train <- sum(colnames(newdata) %in%
                      names(object$wa.optima))
    deshrink <- object$deshrink
    deshrink.fun <- switch(deshrink,
                           inverse = inv.deshrink,
                           classical = class.deshrink,
                           expanded = expand.deshrink,
                           none = no.deshrink)
    X <- object$orig.x
    ENV <- object$orig.env
    ## tolerance options from model
    O <- object$options.tol
    ## Doing CV?
    if(identical(CV, "none")) {
        want <- names(object$wa.optima) %in%
        colnames(newdata)
        want <- names(object$wa.optima)[want]
        if(object$tol.dw) {
            pred <- WATpred(newdata[,want], object$wa.optima[want],
                            object$model.tol[want])
        } else {
            pred <- WApred(newdata[,want], object$wa.optima[want])
        }
        pred <- deshrink.pred(pred, coef(object))
    } else {
        ## CV wanted
        if(identical(CV, "LOO")) {
            loo.pred <- matrix(0, ncol = n.train,
                               nrow = n.fossil)
            mod.pred <- length(n.train)
            useN2 <- object$options.tol$useN2
            want <- names(object$wa.optima) %in% colnames(newdata)
            want <- names(object$wa.optima)[want]
            for(i in seq_len(n.train)) {
                if(verbose && ((i %% 10) == 0)) {
                    cat(paste("Leave one out sample", i, "\n"))
                    flush.console()
                }
                wa.optima <- w.avg(X[-i,], ENV[-i])
                tol <- w.tol(X[-i, ], ENV[-i], wa.optima,
                             useN2 = useN2)
                ## fix up problematic tolerances
                tol <- fixUpTol(tol, O$na.tol, O$small.tol,
                                O$min.tol, O$f, ENV[-i])
                ## CV for the training set
                if(object$tol.dw) {
                    wa.env <- WATpred(X[-i,], wa.optima, tol)
                    mod.pred[i] <- WATpred(X[i,,drop=FALSE],
                                           wa.optima, tol)
                } else {
                    wa.env <- WApred(X[-i,], wa.optima)
                    mod.pred[i] <- WApred(X[i,,drop=FALSE],
                                          wa.optima)
                }
                deshrink.mod <- deshrink.fun(ENV[-i], wa.env)
                wa.env <- deshrink.mod$env
                coefs <- coef(deshrink.mod)
                ## LOO model predictions
                mod.pred[i] <- deshrink.pred(mod.pred[i], coefs)
                ## newdata predictions
                pred <- if(object$tol.dw) {
                    WATpred(newdata[,want], wa.optima[want],
                            tol[want])
                } else {
                    WApred(newdata[,want], wa.optima[want])
                }
                loo.pred[,i] <- deshrink.pred(pred, coefs)
            }
            ## average the LOO predictions
            pred <- rowMeans(loo.pred)
        } else if(identical(CV, "bootstrap")) {
            boot.pred <- matrix(0, ncol = n.boot, nrow = n.fossil)
            oob.pred <- matrix(NA, ncol = n.boot, nrow = n.train)
            for(i in seq_len(n.boot)) {
                if(verbose && ((i %% 100) == 0)) {
                    cat(paste("Bootstrap sample", i, "\n"))
                    flush.console()
                }
                ## bootstrap sample
                sel <- .Internal(sample(n.train, n.train,
                                        TRUE, NULL))
                wa.optima <- w.avg(X[sel,], ENV[sel])
                ## do the model bits
                ones <- rep(1, length = length(wa.optima))
                miss <- is.na(wa.optima)
                ones[miss] <- 0
                wa.optima[miss] <- 0
                rowsum <- X[sel,] %*% ones
                wa.env <- (X[sel,] %*% wa.optima) / rowsum
                deshrink.mod <- deshrink.fun(ENV[sel], wa.env)
                wa.env <- deshrink.mod$env
                coefs <- coef(deshrink.mod) #$coef
                ## if we want sample specific errors or
                ## model performance stats
                rowsum <- X[-sel,] %*% ones
                pred <- (X[-sel,] %*% wa.optima) / rowsum
                oob.pred[-sel,i] <- deshrink.pred(pred, coefs)
                ## do the prediction step
                want <- names(wa.optima) %in% colnames(newdata)
                want <- names(wa.optima)[want]
                ones <- rep(1, length = length(want))
                miss <- miss[want]
                ones[miss] <- 0
                rowsum <- newdata[,want] %*% ones
                pred <- (newdata[,want] %*% wa.optima[want]) /
                    rowsum
                boot.pred[,i] <- deshrink.pred(pred, coefs)
            }
            pred <- rowMeans(boot.pred)
        } else if (identical(CV, "nfold")) {
            boot.pred <- matrix(0, ncol = n.boot, nrow = n.fossil)
            oob.pred <- matrix(NA, ncol = n.boot, nrow = n.train)
            ind <- rep(1:nfold, length = n.train)
            for(i in seq_len(n.boot)) {
                if(verbose && ((i %% 100) == 0)) {
                    cat(paste("n-fold sample", i, "\n"))
                    flush.console()
                }
                ## n-fold sample
                pind <- sample(ind)
                for (k in seq_len(nfold)) {
                    sel <- pind != k
                    wa.optima <- w.avg(X[sel,], ENV[sel])
                    ## do the model bits
                    ones <- rep(1, length = length(wa.optima))
                    miss <- is.na(wa.optima)
                    ones[miss] <- 0
                    wa.optima[miss] <- 0
                    rowsum <- X[sel,] %*% ones
                    wa.env <- (X[sel,] %*% wa.optima) / rowsum
                    deshrink.mod <- deshrink.fun(ENV[sel], wa.env)
                    wa.env <- deshrink.mod$env
                    coefs <- coef(deshrink.mod) #$coef
                    ## if we want sample specific errors or
                    ## model performance stats
                    rowsum <- X[!sel,] %*% ones
                    pred <- (X[!sel,] %*% wa.optima) / rowsum
                    oob.pred[!sel,i] <- deshrink.pred(pred, coefs)
                    ## do the prediction step
                    want <- names(wa.optima) %in% colnames(newdata)
                    want <- names(wa.optima)[want]
                    ones <- rep(1, length = length(want))
                    miss <- miss[want]
                    ones[miss] <- 0
                    rowsum <- newdata[,want] %*% ones
                    pred <- (newdata[,want] %*% wa.optima[want]) /
                        rowsum
                    boot.pred[,i] <- deshrink.pred(pred, coefs)
                }
            }
            pred <- rowMeans(boot.pred)
        }
    }
    .call <- match.call()
    .call[[1]] <- as.name("predict")
    names(pred) <- rownames(newdata)
    retval <- list(pred = list(pred = pred, rmsep = NULL),
                   performance = NULL, model.pred = NULL)
    if(identical(CV, "none")) {
        retval$performance <- with(object,
                                   list(r.squared = r.squared,
                                        avg.bias = avg.bias,
                                        max.bias = max.bias,
                                        rmsep = rmse))
        retval$model.pred <- list(pred = fitted(object))
    } else if(identical(CV, "LOO")) {
        mod.r.squared <- cor(mod.pred, ENV)^2
        mod.resid <- mod.pred - ENV
        mod.avg.bias <- mean(mod.resid)
        mod.max.bias <- maxBias(mod.resid, ENV)
        mod.rmsep <- sqrt(mean(mod.resid^2))
        retval$performance <- list(r.squared = mod.r.squared,
                                   avg.bias = mod.avg.bias,
                                   max.bias = mod.max.bias,
                                   rmsep = mod.rmsep)
        names(mod.pred) <- names(mod.resid) <- rownames(X)
        retval$model.pred <- list(pred = mod.pred,
                                  resid = mod.resid)
    } else {
        mod.pred <- rowMeans(oob.pred, na.rm = TRUE)
        mod.resid <- mod.pred - ENV
        s1 <- apply(oob.pred, 1, sd, na.rm = TRUE)
        s2 <- sqrt(rowMeans((oob.pred - ENV)^2, na.rm = TRUE))
        mod.s1 <- sqrt(mean(s1^2))
        mod.s2 <- sqrt(mean(mod.resid^2))
        samp.rmsep <- sqrt(s1^2 + mod.s2^2)
        mod.rmsep <- sqrt(mod.s1^2 + mod.s2^2)
        mod.r.squared <- cor(mod.pred, ENV)^2
        mod.avg.bias <- mean(mod.resid)
        mod.max.bias <- maxBias(mod.resid, ENV)
        retval$performance <- list(r.squared = mod.r.squared,
                                   avg.bias = mod.avg.bias,
                                   max.bias = mod.max.bias,
                                   rmsep = mod.rmsep)
        names(mod.pred) <- names(mod.resid) <-
            names(samp.rmsep) <- rownames(X)
        retval$model.pred <- list(pred = mod.pred,
                                  resid = mod.resid,
                                  rmsep = samp.rmsep)
        test.s1 <- apply(boot.pred, 1, sd)
        test.rmsep <- sqrt(test.s1^2 + mod.s2^2)
        names(test.rmsep) <- names(retval$pred$pred)
        retval$pred$pred <- pred
        retval$pred$rmsep <- test.rmsep
    }
    retval$call = .call
    if (identical(CV, "nfold"))
        CV <- paste(nfold, "fold", sep="-")
    retval$CV.method <- CV
    retval$deshrink <- deshrink
    retval$tol.dw <- object$tol.dw
    class(retval) <- "predict.wa"
    retval
}

WApred <- function(X, optima) {
    ones <- rep.int(1, length(optima))
    miss <- is.na(optima)
    ones[miss] <- 0
    optima[miss] <- 0
    rsum <- X %*% ones
    ((X %*% optima) / rsum)
}

WATpred <- function(X, optima, tol) {
    ones <- rep.int(1, length(optima))
    miss <- is.na(optima)
    ones[miss] <- 0
    optima[miss] <- 0
    tol[miss] <- 1
    tol2 <- tol^2
    #tmp <- sweep(X, 2, optima, "*", check.margin = FALSE)
    ##tmp <- RowSums(t(t(X) * optima / tol2))
    #tmp <- RowSums(sweep(tmp, 2, tol2, "/",
    #                     check.margin = FALSE))
    #tmp / RowSums(sweep(X, 2, tol2, "/",
    #                    check.margin = FALSE)[,!miss, drop = FALSE])
    #tmp / (sweep(X, 2, tol2, "/",check.margin = FALSE) %*% ones)
    ##tmp / (t(t(X) / tol2) %*% ones)
    RowSums(t(t(X) * optima / tol2)) / (t(t(X) / tol2) %*% ones)
}
