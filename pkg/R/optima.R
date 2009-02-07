## User function to compute weighted average optima for taxa

`optima` <- function(x, ...)
    UseMethod("optima")

`optima.default` <- function(x, env, ...) {
    x <- data.matrix(x)
    opt <- colSums(x * env) / colSums(x)
    names(opt) <- colnames(x)
    class(opt) <- "optima"
    attr(opt, "env") <- deparse(substitute(env))
    opt
}

`print.optima` <- function(x, ...) {
    cat("\n")
    msg <- paste("Weighted Average Optima For:", attr(x, "env"))
    writeLines(strwrap(msg, prefix = "\t"),
               sep = "\n\n")
    attr(x, "env") <- NULL
    print(unclass(x), ...)
}
