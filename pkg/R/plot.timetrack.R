`plot.timetrack` <- function(x, choices = 1:2, order,
                             ptype = c("l", "p", "o", "b"),
                             pch = c(1,2),
                             col = c("black","red"),
                             lty = "solid", lwd = 1,
                             ...) {
    ptype <- match.arg(ptype)
    plt <- plot(x$ord, choices = choices, scaling = x$scaling,
                type = "p", display = "sites", ...,
                pch = pch[1], col = col[1])
    pass <- fitted(x, type = "passive", choices = choices)
    if(!missing(order)) {
        if(length(order) != NROW(pass))
            stop("'length(order)' not equal to number of passive samples.")
        pass[order, ]
    }
    if(ptype %in% c("l", "o", "b")) {
        lines(pass, pch = pch[2], col = col[2],
              lty = lty, lwd = lwd, type = ptype, ...)
    } else {
        points(pass, pch = pch[2], col = col[2], ...)
    }
    invisible()
}
