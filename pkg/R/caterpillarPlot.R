`caterpillarPlot` <- function(x, env, useN2 = TRUE, decreasing = TRUE,
                              mult = 1, labels, xlab = NULL,
                              pch = 21, bg = "white", col = "black", lcol = col,
                              ...) {
    ## compute the optima
    opt <- optima(x = x, env = env)
    ## and tolerances
    tol <- tolerance(x = x, env = env, useN2 = useN2)

    ## reorder
    ord <- order(opt, decreasing = decreasing)
    opt <- opt[ord]
    tol <- tol[ord]

    ## par
    op <- par(yaxs = "i")
    on.exit(par(op))

    ## number of species
    nspp <- ncol(x)
    yvals <- seq_len(nspp)

    ## labels == spp names
    if(missing(labels)) {
        labels <- names(opt)
    }
    linch <- if (!is.null(labels))
        max(strwidth(labels, "inch"), na.rm = TRUE)
    nmai <- par("mai")
    nmai[2L] <- nmai[4L] + linch + 0.1
    par(mai = nmai)

    ## ylab
    if(is.null(xlab))
        xlab <- deparse(substitute(env))

    ## tolerance range
    upr <- opt + (tol * mult)
    lwr <- opt - (tol * mult)

    ## do the plot
    plot(c(lwr, upr), rep.int(yvals, 2), type = "n", axes = FALSE,
         ylab = "", xlab = xlab, ylim = range(0, yvals + 1), ...)
    abline(h = yvals, lty = 1, lwd = 0.5, col = "lightgray")
    segments(lwr, yvals, upr, yvals, col = lcol, ...)
    points(opt, yvals, pch = pch, bg = bg, col = col, ...)
    axis(side = 1, ...)
    axis(side = 2, labels = labels, at = yvals, las = 1, ...)
    box()
    out <- data.frame(Optima = opt, Tolerance = tol)
    invisible(out)
}
