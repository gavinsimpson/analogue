`n2` <- function(x, ...) {
    UseMethod("n2")
}

`n2.default` <- function(x, which = c("species", "sites"), ...) {
    which <- match.arg(which)
    if (isTRUE(all.equal(which, "species"))) {
        MARGIN <- 2
        sumFun <- colSums
    } else {
        MARGIN <- 1
        sumFun <- rowSums
    }
    pj <- sumFun(x)
    pj <- sumFun(sweep(x, MARGIN, pj, "/")^2, na.rm = TRUE)
    N2 <- 1 / pj
    N2
}
