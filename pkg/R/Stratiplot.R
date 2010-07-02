## Stratigraphic plots using lattice
`Stratiplot` <- function(x, ...)
    UseMethod("Stratiplot")

`Stratiplot.default` <- function(x, y,
                                 type = "l",
                                 ylab = NULL,
                                 xlab = "",
                                 pages = 1,
                                 rev = TRUE,
                                 ylim,
                                 sort = c("none", "wa", "var"),
                                 svar = NULL,
                                 rev.sort = FALSE,
                                 strip = FALSE,
                                 ...) {
    ## inline function for custom axis
    axis.VarLabs <- function(side, ...) {
        if(isTRUE(all.equal(side, "top"))) {
            M <- function(lims) min(lims) + (diff(lims) / 2)
            xlim <- current.panel.limits()$xlim
            panel.axis(side = side, outside = TRUE, at = 0, #M(xlim),
                       tck = 1, line.col = "black",
                       text.col = "black",
                       labels = levels(sx$ind)[which.packet()],
                       rot = 60)
        } else {
            axis.default(side = side, ...)
        }
    }

    ## process 'type'
    TYPE <- c("l","h","g","smooth","b","o", "poly","p")
    if(any(!(type %in% TYPE)))
        stop("Invalid 'type' specified")
    n.vars <- NCOL(x)
    ## ylabel
    if(is.null(ylab))
        ylab <- deparse(substitute(y))
    ## do we need to sort variables
    if(missing(sort))
        sort <- "none"
    sort <- match.arg(sort)
    if((check.var <- (missing(svar) || is.null(svar))))
        svar <- y
    ord <- seq_len(n.vars)
    if(sort == "wa") {
        ## sort by
        opt <- optima(x, svar)
        ord <- order(opt)
    } else if(sort == "var") {
        if(check.var) {
            warning("With 'sort = \"var\"', 'svar' not supplied.\nNo sorting applied.")
        } else {
            ord <- order(svar)
        }
    }
    if(rev.sort)
        ord <- rev(ord)
    x <- x[, ord]
    sx <- stack(x)
    sx$ind <- factor(sx$ind, levels = colnames(x))
    ## check length of y
    if(!isTRUE(all.equal((leny <- length(y)), (nr <- nrow(sx))))) {
        ## if length(y) == nrow(sx)/n.vars, then expand
        if(isTRUE(all.equal(leny, nr / n.vars)))
            y <- rep(y, n.vars)
        else
            stop("Ambiguous 'length(y)';\nmust be equal to 'nrow(x)' or\n'nrow(x) / number of species'.")
    }
    ## plot parameters
    maxy <- max(y)
    miny <- min(y)
    ## add padYlim * range as per base graphics
    padY <- 0.01
    if(missing(ylim)) {
        diffy <- padY * (maxy - miny)
        ylim <- c(miny - diffy, maxy + diffy)
    } else {
        minLim <- min(ylim)
        maxLim <- max(ylim)
        ## add padY * range as per base graphics
        diffy <- abs(diff(minLim, maxLim))
        ylim <- if(minLim > maxLim)
            c(minLim + (padY * diffy), maxLim - (padY * diffy))
        else
            c(minLim - (padY * diffy), maxLim + (padY * diffy))
    }
    if(rev)
        ylim <- rev(ylim)
    max.abun <- sapply(x, function(x) round(max(x), 1))
    xlimits <- lapply(max.abun * 1.05, function(x) c(0, x))
    scales <- list(cex = 0.75, tck = 0.75,
                   y = list(axs = "r", limits = ylim),
                   x = list(axs = "r", rot = 45, relation = "free"))
    par.strip.text <- list(cex = 0.75)
    str.max <- 1
    if(!isTRUE(strip)) {
        gp <- gpar()
        convWidth <- function(x, gp) {
            convertWidth(grobWidth(textGrob(x, gp = gp)), "lines",
                         valueOnly = TRUE)
        }
        str.max <- max(sapply(levels(sx$ind), convWidth, gp))
        str.max <- ceiling(str.max) + 6
    }
    ## plotting
    xyplot(y ~ values | ind,
           data = sx,
           type = type,
           ylab = ylab, xlab = xlab,
           strip.left = FALSE, strip = strip,
           par.strip.text = par.strip.text,
           scales = scales,
           xlim = xlimits,
           ylim = ylim,
           panel = "panel.Stratiplot",
           layout = c(n.vars, 1, pages),
           par.settings = list(layout.widths = list(panel = max.abun),
           layout.heights = list(top.padding = str.max)),
           axis = if(isTRUE(strip)) {axis.default} else {axis.VarLabs},
           ...)
}
