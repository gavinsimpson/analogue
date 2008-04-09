`print.performance` <- function(x, ...) {
    CV.method <- attr(x, "CV.method")
    perf.names <- names(x)
    attributes(x) <- NULL
    names(x) <- perf.names
    print.default(round(x, 4), ...)
    invisible(x)
}
