`plot.timetrack` <- function(x, choices = 1:2,
                             pch = c(1,2),
                             col = c("black","red"),
                             ...) {
    plt <- plot(x$ord, choices = choices, scaling = x$scaling,
                type = "p", display = "sites", ...,
                pch = pch[1], col = col[1])
    ##points(fitted(x)[, choices])
    lines(fitted(x)[, choices], pch = pch[2], col = col[2])
}
