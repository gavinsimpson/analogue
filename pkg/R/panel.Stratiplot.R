`panel.Stratiplot` <- function(x, y, type = "l",
                               col,
                               pch = plot.symbol$pch,
                               cex = plot.symbol$cex,
                               col.line = plot.line$col,
                               col.symbol = plot.symbol$col,
                               col.refline = ref.line$col,
                               col.smooth = "red",
                               col.poly = plot.line$col,
                               lty = plot.line$lty,
                               lwd = plot.line$lwd,
                               lty.smooth = plot.line$lty,
                               lwd.smooth = 2,
                               fill = plot.symbol$fill,
                               ...) {
    if (all(is.na(x) | is.na(y)))
        return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    plot.symbol <- trellis.par.get("plot.symbol")
    plot.line <- trellis.par.get("plot.line")
    ref.line <- trellis.par.get("reference.line")
    if (!missing(col)) {
        if (missing(col.line))
            col.line <- col
        if (missing(col.symbol))
            col.symbol <- col
    }
    panel.refline(v = 0, col.line = ref.line, ...)
    if ("o" %in% type || "b" %in% type)
        type <- c(type, "p", "l")
    if ("g" %in% type)
        panel.grid(h = -1, v = -1)
    if("l" %in% type)
        panel.lines(x = x, y = y, col = col.line,
                    lty = lty, lwd = lwd, ...)
    if("p" %in% type)
        panel.points(x = x, y = y, cex = cex, fill = fill,
                     col = col.symbol, pch = pch, ...)
    if("h" %in% type)
        panel.segments(x0 = 0, y0 = y, x1 = x, y1 = y,
                       col = col.line, lty = lty, lwd = lwd, ...)
    if("poly" %in% type)
        panel.polygon(x = c(0, x, 0), y = c(y[1], y, y[length(y)]),
                      border = col.poly, col = col.poly, ...)
    if("smooth" %in% type)
        panel.Loess(x, y, col = col.smooth, lwd = lwd.smooth,
                    lty = lty.smooth, ...)
}
