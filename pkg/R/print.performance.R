`print.performance` <- function(x, digits = min(getOption("digits"), 4),
                                ...) {
    CV.method <- attr(x, "CV.method")
    if(inherits(x, "data.frame")) {
        print.data.frame(x, digits = digits, ...)
    } else {
        perf.names <- names(x)
        attributes(x) <- NULL
        names(x) <- perf.names
        print.default(x, digits = digits, ...) # x was round(x, 4)
    }
    invisible(x)
}
