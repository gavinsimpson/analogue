## User function to compute weighted average tolerances for taxa

`tolerance` <- function(x, ...)
    UseMethod("tolerance")

`tolerance.default` <- function(x, env, useN2 = TRUE, ...) {
    x <- data.matrix(x)
    opt <- optima(x, env, ...)
    tol <- sqrt(colSums(x * outer(env, opt, "-")^2) / colSums(x))
    if(useN2) {
        N2 <- sppN2(x)
        tol <- tol / sqrt(1 - (1 / N2))
    }
    names(tol) <- colnames(x)
    tol
}
