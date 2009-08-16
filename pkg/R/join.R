###########################################################################
##                                                                       ##
## join() - Function to merge any number of data frames                  ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.5-0                                                 ##
## Last modified : 16-Aug-2009                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## ...                  - the data frames to be merged                   ##
## verbose = TRUE       - prints a summary of the number of rows/cols in ##
##                        x, y and the merged data frame                 ##
## na.replace = TRUE    - replaces NA's with zeroes (0) in the merged    ##
##                        data frame                                     ##
##                                                                       ##
## HISTORY:                                                              ##
## 17-Apr-2006 - GLS - 0.1-1 * Function created                          ##
## 31-May-2006 - GLS - 0.1-2 * Original function failed if there were    ##
##                             matching rows.                            ##
##                           * Removed call to merge()                   ##
##                           * Added solution provided by Sundar Dorai-  ##
##                             Raj, as indicated.                        ##
##                           * Added some error checking, assume '...'   ##
##                             are all data frames.                      ##
##                           * Updated the verbose sections to match     ##
##                           * Updated documentation                     ##
## 05-Jul-2006 - GLS - 0.1-3 * join() was dropping the rownames of the   ##
##                             joined objects. FIXED                     ##
## 23-Jul-2007 - GLS - 0.2-0 * More user-friendly; will unsplit joined   ##
##                             datasets if split == TRUE                 ##
## 13-Oct-2007 - GLS - 0.3-0 * join() now merges factors correctly       ##
## 15-Oct-2007 - GLS - 0.4-0 * join() now returns a classed object       ##
## 27-Apr-2008 - GLS - 0.4-0 * join() now checks for inheritance from    ##
##                             class data.frame, not that it is that     ##
##                             class. This allows join to work with the  ##
##                             results of join(..., split = FALSE        ##
## 16-Aug-2009 - GLS - 0.5-0 * now has left outer as well as outer join  ##
##                           * replacement value can now be specified    ##
##                                                                       ##
###########################################################################
join <- function(..., verbose = FALSE, na.replace = TRUE, split = TRUE,
                  value = 0, type = c("outer","left"))
{
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
        joined <- do.call(rbind, X)
        colnames(joined) <- cn
        return(joined)
    }
    leftJoin <- function(X) {
        cn <- unique(unlist(lapply(X, colnames)[[1]]))
        ## if more than 2 df in X, merge all bar first
        if(length(X) > 2)
            dfs <- outerJoin(X[-1])
        else
            dfs <- X[[2]]
        ## matched column names
        mcn <- match(colnames(dfs), cn)
        mcn <- mcn[!is.na(mcn)]
        joined <- matrix(NA, ncol = dims[1,2], nrow = sum(dims[,1]))
        joined[1:dims[1,1], ] <- data.matrix(X[[1]])
        joined[(dims[1,1]+1):NROW(joined), mcn] <- data.matrix(dfs[, mcn])
        colnames(joined) <- cn
        return(joined)
    }
    x <- list(...)
    if(any(!sapply(x, inherits, "data.frame")))
        stop("\nall objects to be merged must be data frames.")
    dims <- do.call(rbind, lapply(x, dim))
    n.joined <- nrow(dims)
    if(missing(type))
        type <- "outer"
    type <- match.arg(type)
    if(type == "outer") {
        joined <- outerJoin(x)
    } else if(type == "left") {
        joined <- leftJoin(x)
    }
    if(na.replace) {
        joined[is.na(joined)] <- value
    }
    rn <- lapply(x, rownames)
    if(verbose) {
        stats <- rbind(dims, dim(joined))
        rownames(stats) <- c(paste("Data set ", c(1:n.joined), ":", sep = ""),
                             "Merged:")
        colnames(stats) <- c("Rows", "Cols")
        cat("\nSummary:\n\n")
        printCoefmat(stats, digits = max(3, getOption("digits") - 3),
                     na.print = "")
        cat("\n")
    }
    if(split) {
        retval <- vector(mode = "list", length = n.joined)
        ends<- cumsum(dims[,1])
        start <- c(1, ends[-n.joined] + 1)
        for(i in 1:n.joined) {
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
