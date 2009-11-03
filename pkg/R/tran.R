`tran` <- function(x, ...) {
    UseMethod("tran")
}

`tran.default` <- function(x, method, a = 1, b = 0, p = 2, base = exp(1),
                   na.rm = FALSE, na.value = 0, ...) {
    wasDF <- is.data.frame(x)
    dim.nams <- dimnames(x)
    x <- data.matrix(x)
    METHOD <- c("sqrt", "cubert", "rootroot", "log", "reciprocal", "freq",
                "center", "standardize", "range", "percent", "proportion",
                "pa","missing", "hellinger", "chi.square", "wisconsin",
                "pcent2prop", "prop2pcent", "logRatio", "power",
                "rowCenter")
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
                    cubert = sign(x) * exp(log(abs(x)) / 3), #x^(1/3),
                    rootroot = sign(x) * exp(log(abs(x)) / 3), #x^(1/4),
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
                    prop2pcent = x * 100,
                    logRatio = {x <- sweep(x, 2, a, "*")
                                x <- sweep(x, 2, b, "+")
                                x <- log(x, base = base)
                                x - rowMeans(x)},
                    power = x^p,
                    rowCenter = x - rowMeans(x)
                    )
    }
    if(wasDF)
        x <- as.data.frame(x)
    dimnames(x) <- dim.nams
    attr(x, "tran") <- method
    return(x)
}

`tran.formula` <- function(formula, data = NULL,
                           subset = NULL,
                           na.action = na.pass, ...) {
    mf <- match.call()
    mf[[1]] <- as.name("model.frame")
    mt <- terms(formula, data = data, simplify = TRUE)
    mf[[2]] <- formula(mt, data = data)
    mf$na.action <- substitute(na.action)
    dots <- list(...)
    mf[names(dots)] <- NULL
    mf <- eval(mf, parent.frame())
    tran.default(mf, ...)
}
