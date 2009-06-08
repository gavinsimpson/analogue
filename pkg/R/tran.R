`tran` <- function(x, method, a = 1, b = 0, base = exp(1),
                   na.rm = FALSE, na.value = 0, ...) {
    wasDF <- is.data.frame(x)
    dim.nams <- dimnames(x)
    x <- data.matrix(x)
    METHOD <- c("sqrt", "cubert", "log", "reciprocal", "freq", "center",
                "standardize", "range", "percent", "proportion", "pa",
                "missing", "hellinger", "chi.square", "wisconsin",
                "pcent2prop", "prop2pcent")
    method <- match.arg(method, METHOD)
    if(method %in% c("freq", "standardize","range","pa","hellinger",
                     "chi.square","wisconsin")) {
        if(isTRUE(all.equal(method, "wisconsin")))
            x <- wisconsin(x)
        else
            x <- decostand(x, method = method, na.rm = na.rm, ...)
        attr(x, "decostand") <- NULL
    } else {
        x <- switch(method,
                    sqrt = sqrt(x),
                    cubert = x^(1/3),
                    log = {x <- sweep(x, 2, a, "*")
                           x <- sweep(x, 2, b, "+")
                           log(x, base = base)} ,
                    reciprocal = 1 / x,
                    center = scale(x, scale = FALSE, center = TRUE),
                    percent = sweep(x, 1, rowSums(x), "/") * 100,
                    proportion = sweep(x, 1, rowSums(x), "/"),
                    missing = apply(x, 2,
                    function(x) {x[is.na(x)] <- na.value
                                 x}),
                    pcent2prop = x / 100,
                    prop2pcent = x * 100
                    )
    }
    if(wasDF)
        x <- as.data.frame(x)
    dimnames(x) <- dim.nams
    attr(x, "tran") <- method
    return(x)
}
