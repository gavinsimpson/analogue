## compute the standard error of MAT reconstructed values
## following the ideas of ter Braak 1995 to use the weighted
## variance of the k-closest analogues

`stdError` <- function(object, ...) {
    UseMethod("stdError")
}

`stdError.mat` <- function(object, k, ...) {
    getOrd <- function(dis) {
        nas <- is.na(dis)
        order(dis[!nas])
    }
    getWts <- function(i, dis, ords, k.seq) {
        nas <- is.na(dis[,i])
        dis[!nas, i][ords[,i]][k.seq]
    }
    getEnv <- function(i, dis, ords, k.seq, y){
        nas <- is.na(dis[,i])
        y[!nas][ords[,i]][k.seq]
    }
    if(missing(k)) {
        k <- getK(object, ...)
        auto <- object$auto
    } else {
        auto <- FALSE
    }
    ## create k sequence
    k.seq <- seq_len(k)
    ## ordering of objects in terms of dissim
    ords <- apply(object$Dij, 2, getOrd)
    SEQ <- seq_len(ncol(ords))
    ## weights = 1/Dij
    wts <- 1 / sapply(SEQ, getWts, object$Dij, ords, k.seq)
    ## produce matrix of Env data for each site
    env <- sapply(SEQ, getEnv, object$Dij, ords, k.seq,
                  object$orig.y)
    ## mean of env of k closest analogues
    ybar <- colMeans(env)
    wtdSD <- sqrt(colSums(wts * sweep(env, 2, ybar, "-")^2) /
                  colSums(wts))
    names(wtdSD) <- names(object$orig.y)
    class(wtdSD) <- "stdError"
    return(wtdSD)
}

`stdError.predict.mat` <- function(object, k, ...) {
    getOrd <- function(dis) {
        nas <- is.na(dis)
        order(dis[!nas])
    }
    getWts <- function(i, dis, ords, k.seq) {
        nas <- is.na(dis[,i])
        dis[!nas, i][ords[,i]][k.seq]
    }
    getEnv <- function(i, dis, ords, k.seq, y){
        nas <- is.na(dis[,i])
        y[!nas][ords[,i]][k.seq]
    }
    if(missing(k)) {
        k <- getK(object, ...)
        auto <- object$auto
    } else {
        auto <- FALSE
    }
    ## create k sequence
    k.seq <- seq_len(k)
    ## ordering of objects in terms of dissim
    ords <- apply(object$Dij, 2, getOrd)
    SEQ <- seq_len(ncol(ords))
    ## weights = 1/Dij
    wts <- 1 / sapply(SEQ, getWts, object$Dij, ords, k.seq)
    ## produce matrix of Env data for each site
    env <- sapply(SEQ, getEnv, object$Dij, ords, k.seq,
                  object$observed)
    ## mean of env of k closest analogues
    ybar <- colMeans(env)
    wtdSD <- sqrt(colSums(wts * sweep(env, 2, ybar, "-")^2) /
                  colSums(wts))
    names(wtdSD) <- colnames(object$predictions$model$predicted)
    class(wtdSD) <- "stdError"
    attr(wtdSD, "k") <- k
    attr(wtdSD, "auto") <- object$auto
    return(wtdSD)
}
