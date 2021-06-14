`join` <- function(..., verbose = FALSE, na.replace = TRUE,
                   split = TRUE, value = 0,
                   type = c("outer", "left", "inner")) {
    outerJoin <- function(X) {
        ## From code provided by Sundar Dorai-Raj in R-Help posting:
        ## http://article.gmane.org/gmane.comp.lang.r.general/63042/
        cn <- unique(unlist(lapply(X, colnames)))
        for(i in seq(along = X)) {
            if(any(m <- !cn %in% colnames(X[[i]]))) {
                na <- matrix(NA, nrow(X[[i]]), sum(m))
                dimnames(na) <- list(rownames(X[[i]]), cn[m])
                X[[i]] <- cbind(X[[i]], na)
            }
        }
        joined <- do.call("rbind", X)
        colnames(joined) <- cn
        joined
    }
    leftJoin <- function(X, dims) {
        cn <- colnames(X[[1L]])
        ## if more than 2 df in X, merge all bar first
        if(length(X) > 2) {
            dfs <- outerJoin(X[-1])
        } else {
            dfs <- X[[2]]
        }
        ## matched column names
        mcn <- match(colnames(dfs), cn)
        mcn <- mcn[!is.na(mcn)]         # remove NAs
        joined <- matrix(NA, ncol = dims[1L, 2L], nrow = sum(dims[, 1L]))
        colnames(joined) <- cn
        joined[seq_len(dims[1L, 1L]), ] <- data.matrix(X[[1L]])
        joined[(dims[1L, 1L]+1L):NROW(joined), mcn] <- data.matrix(dfs[, cn[mcn]])
        joined
    }
    innerJoin <- function(X) {
        cn <- lapply(X, colnames)
        cn <- Reduce(intersect, cn)
        ##joined <- matrix(NA, ncol = length(cn), nrow = sum(dims[,1]))
        joined <- vector(length = length(X), mode = "list")
        for(i in seq_along(joined)) {
            joined[[i]] <- data.matrix(X[[i]][, cn])
        }
        joined <- do.call("rbind", joined)
        colnames(joined) <- cn
        joined
    }
    x <- list(...)
    if(any(!sapply(x, inherits, "data.frame", USE.NAMES = FALSE)))
        stop("\nall objects to be merged must be data frames.")
    dims <- t(vapply(x, FUN = dim, FUN.VALUE = integer(2)))
    n.joined <- nrow(dims)
    if(missing(type)) {
        type <- "outer"
    }
    type <- match.arg(type)
    joined <- switch(type,
                     outer = outerJoin(x),
                     left  = leftJoin(x, dims = dims),
                     inner = innerJoin(x))
    if(na.replace) {
        joined[is.na(joined)] <- value
    }
    rn <- lapply(x, rownames)
    if(verbose) {
        stats <- rbind(dims, dim(joined))
        rownames(stats) <- c(paste("Data set ", seq_len(n.joined), ":", sep = ""),
                             "Merged:")
        colnames(stats) <- c("Rows", "Cols")
        cat("\nSummary:\n\n")
        printCoefmat(stats, digits = 3,
                     na.print = "")
        cat("\n")
    }
    if(split) {
        retval <- vector(mode = "list", length = n.joined)
        ends<- cumsum(dims[,1])
        start <- c(1, ends[-n.joined] + 1)
        for(i in seq_len(n.joined)) {
            retval[[i]] <- as.data.frame(joined[start[i]:ends[i], ])
            rownames(retval[[i]]) <- rn[[i]]
        }
        names(retval) <- as.character(match.call())[2:(n.joined+1)]
        class(retval) <- "join"
    } else {
        retval <- as.data.frame(joined, row.names = rownames(joined))
        class(retval) <- c("join", class(retval))
    }
    attr(retval, "type") <- type
    retval
}
