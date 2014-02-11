`plot.timetrack` <- function(x, choices = 1:2,
                             display = c("wa","lc"),
                             order,
                             ptype = c("l", "p", "o", "b"),
                             pch = c(1,2),
                             col = c("black","red"),
                             lty = "solid", lwd = 1,
                             ...) {
    ptype <- match.arg(ptype)
    display <- match.arg(display)
    scrs <- scores(x$ord, choices = choices, scaling = x$scaling,
                   display = display, ...)
    pass <- fitted(x, type = "passive", choices = choices)
    xlim <- range(scrs[,1], pass[,1])
    ylim <- range(scrs[,2], pass[,2])
    plt <- plot(x$ord, choices = choices, scaling = x$scaling,
                type = "n", display = display, ...,
                ylim = ylim, xlim = xlim)
    points(scrs, pch = pch[1], col = col[1], ...)
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
