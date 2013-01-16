`fitted.timetrack` <-
    function(object, type = c("passive", "ordination"),
             model = NULL, choices = 1:2, ...)
{
    if(missing(type))
        type <- "passive"
    type <- match.arg(type)
    model <- if(is.null(model)) {
        if(is.null(object$ordination$CCA)) "CA" else "CCA"
    }
    if(isTRUE(all.equal(type, "passive"))) {
        fit <- fitted(unclass(object), ...)[, choices, drop = FALSE]
    } else {
        fit <- fitted(object$ordination, model = model,
                      ...)[, choices, drop = FALSE]
    }
    fit
}
