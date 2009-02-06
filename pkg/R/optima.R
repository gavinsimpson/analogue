## User function to compute weighted average optima for taxa

`optima` <- function(x, ...)
    UseMethod("optima")

`optima.default` <- function(x, env, ...) {
    x <- data.matrix(x)
    opt <- colSums(x * env) / colSums(x)
    names(opt) <- colnames(x)
    opt
}
