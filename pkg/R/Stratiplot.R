## Stratigraphic plots using lattice
`Stratiplot` <- function(x, ...)
    UseMethod("Stratiplot")

`Stratiplot.default` <- function(x, y,
                                 type = "l",
                                 ylab = NULL,
                                 xlab = "",
                                 pages = 1,
                                 rev = TRUE,
                                 sort = c("none", "wa", "var"),
                                 svar = NULL,
                                 rev.sort = FALSE,
                                 ...) {
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
    } else if(check.var) {
        warning("With 'sort = \"var\"', 'svar' not supplied.\nNo sorting applied.")
    } else {
        ord <- order(svar)
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
    ylimits <- c(0 - (0.03*maxy), maxy + (0.03 * maxy))
    if(rev)
        ylimits <- rev(ylimits)
    max.abun <- sapply(x, function(x) round(max(x), 1))
    xlimits <- lapply(max.abun * 1.05, function(x) c(0, x))
    scales <- list(cex = 0.75, tck = 0.75,
                   y = list(axs = "r", limits = ylimits),
                   x = list(axs = "r", rot = 45, relation = "free"))
    par.strip.text <- list(cex = 0.75)
    ## plotting
    xyplot(y ~ values | ind,
           data = sx,
           type = type,
           ylab = ylab, xlab = xlab,
           strip.left = FALSE, strip = TRUE,
           par.strip.text = par.strip.text,
           scales = scales,
           xlim = xlimits,
           prepanel = function(x, y, ...) {
               list(xlim = c(0, max(x)),
                    ylim = rev(c(0, max(y))))
           },
           panel = panel.Stratiplot,
           layout = c(n.vars, 1, pages),
           par.settings = list(layout.widths = list(panel = max.abun)),
           ##index.cond = list(ord),
           ...)
}
