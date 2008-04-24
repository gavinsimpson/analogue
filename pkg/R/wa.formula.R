`wa.formula` <- function(formula, data, subset, na.action,
                         deshrink = c("inverse", "classical", "expanded", "none"),
                         tol.dw = FALSE, model = FALSE, ...) {
    ## set default deshrinking to inverse if no supplied
    if(missing(deshrink))
        deshrink <- "inverse"
    deshrink <- match.arg(deshrink)
    ## the function call
    .call <- match.call()
    ## need to reset due to method dispatch
    .call[[1]] <- as.name("wa")
    ## keep only the arguments which should go into the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## drop the intercept
    attr(attr(mf, "terms"), "intercept") <- 0
    ## 1) allow model.frame to update the terms object before
    ##    saving it.
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    res <- wa.default(x, y, deshrink = deshrink, tol.dw = tol.dw)
    res$na.action <- attr(mf, "na.action")
    res$call <- .call
    if(model) {
        res$terms <- mt
        res$model <- mf
    }
    return(res)
}
