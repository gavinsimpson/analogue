## gradientDist - function for finding distance along gradient
## given object (a vector or something coercible to) that
## represents the gradient position, standardise to locations
## into the range 0, ..., 1
gradientDist <- function(object, ...) {
    UseMethod("gradientDist")
}

gradientDist.default <- function(object, order, na.rm = TRUE, ...) {
    object <- as.vector(object)
    if(missing(order))
        order <- seq_along(object)
    object <- object[order]
    minD <- min(object, na.rm = na.rm)
    k <- if(any(object < 0, na.rm = na.rm)) {
        minD
    } else {
        .Machine$double.eps
    }
    ran <- max(object, na.rm = na.rm)
    ran <- ran - minD
    ran <- pmax(k, ran, na.rm = na.rm)
    object <- object - minD
    object <- object / ran
    class(object) <- "gradientDist"
    object
}

gradientDist.cca <- function(object, na.rm = TRUE, axis = 1L,
                             scaling = 0, ...) {
    if(length(axis) > 1L) {
        axis <- axis[1L]
    }
    scrs <- as.vector(scores(object, choices = axis, scaling = scaling,
                             display = "sites", ...))
    gradientDist.default(scrs, na.rm = na.rm, ...)
}

gradientDist.prcurve <- function(object, na.rm = TRUE, ...) {
    order <- object$tag
    scrs <- object$lambda
    gradientDist.default(scrs, order, na.rm = na.rm, ...)
}
