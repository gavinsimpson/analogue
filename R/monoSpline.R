##' Fits a cubic spline to relate sediment core depth to sediment age, using a monotonicity constraint to maintain stratigrphic ordering in time.
##'
##' @title Monotonic cubic spline age models
##' @param depth numeric; vector of sediment depths
##' @param date numeric; vector of sediment ages/years
##' @param error numeric; vector of errors on the estimated ages
##' @param increasing logical; should the monotonicity be increasing or decreasing (\code{FALSE}, the default). Depending on what you pass as argument \code{date}, you may need to set \code{increasing = TRUE}, say for years before present.
##' @param bs character; the basis type to use. Not sure it makes sense to use other basis types available in \pkg{mgcv}.
##' @param k numeric; the dimension of the basis expansion.
##' @param fx logical; use a fixed degree of freedom spline or use penalized regression to estimate the degree of smoothness.
##' @param ... additional arguments passed to other methods
##' @return An object of class \code{monoSpline} a list with the following components:
##' \item{y}{the observed sample ages or dates as supplied to argument \code{date}}
##' \item{x}{the observed sample depths as supplied to argument \code{depths}}
##' \item{fitted.values}{numeric vector of fitted ages/dates}
##' \item{coefficients}{vector of coefficients for the monotonic spline}
##' \item{residuals}{vector of residuals}
##' \item{sm}{smoothing matrix. See \code{\link[mgcv]{gamObject}}.}
##'
##' @author Gavin L. Simpson
`monoSpline` <- function(depth, date, error = NULL, increasing = FALSE,
                         bs = "cr", k = 10, fx = FALSE, ...) {
    require("mgcv")
    if(is.null(error)) {
        error <- rep(0, length(depth))
        W <- rep(1, length(depth))
    } else {
        W <- 1 / error
        W[1] <- W[2] * error[2]
    }
    df <- data.frame(depth = depth, date = date, error = error)
    mod <- gam(date ~ s(depth, k = k, bs = bs, fx = fx), data = df,
               fit = FALSE)
    sm <- smoothCon(s(depth, k = k, bs = bs, fx = fx), data = df,
                    knots = NULL)[[1]]
    ## Fm are the constraints to enforce monotonicity
    Fm <- mono.con(sm$xp, up = increasing)
    if (!increasing) {
        Xp <- -sm$xp
    } else {
        Xp <- sm$xp
    }
    G <- list(X = sm$X, C = matrix(0,0,0), sp = mod$sp, p = Xp,
              y = df$date, Ain = Fm$A, bin = Fm$b,
              S = sm$S, off = 0, w = W)
    ## fitted parameters
    p <- pcls(G)
    ## fitted values
    pm <- Predict.matrix(sm, data = data.frame(depth = df$depth))
    fit <- pm %*% p
    ## set up return object
    ## this needs all info for a predict() method to work, inc
    ## sm - smoothing matrix
    ## p  - a vector of fitted coefficients; the object returned from pcls
    ##
    ## also return fitted values
    out <- list(y = df$date,
                x = df$depth,
                fitted.values = drop(fit),
                coefficients = p,
                residuals = drop(df$date - fit),
                sm = sm)
    class(out) <- "monoSpline"
    out
}

`predict.monoSpline` <- function(object, newdata, ...) {
    if(missing(newdata))
        yhat <- fitted(object)
    else {
        yhat <- Predict.matrix(object$sm, data = newdata) %*% coef(object)
    }
    yhat
}

`extrap` <- function(obj, to, ...) {
    stopifnot(!missing(to))
    nr <- nrow(obj)
    r <- obj$date[nr] - obj$date[nr-1]
    d <- obj$depth[nr] - obj$depth[nr-1]
    dates <- seq(from = obj$date[nr] + r, to = to, by = r)
    depths <- seq(from = obj$depth[nr] + d, by = d,
                  length.out = length(dates))
    rbind(obj, cbind(depth = depths, date = dates))
}

## Example usage
## dat <- read.csv("/home/gavin/work/projects/greenland/data/dating/sample-dating-data.csv")
## head(dat)

## mod <- with(dat, monoSpline(Depth, Date, Error))

## ## plot observed and fitted
## plot(fitted(mod), mod$y, panel.first = abline(0, 1, col = "grey"))

## ## predict for all levels, here core was in 0.25 slices
## newdat <- data.frame(Depth = seq(0, max(dat$Depth) + 0.25, by = 0.25))
## ## I would probably make these depths be the mid points, but that requires
## ## a little more thought, something like
## newdat2 <- transform(newdat, Depth = Depth - c(0, diff(Depth) / 2))

## ## predict for all core levels
## pred <- transform(newdat2, Year = predict(mod, newdata = newdat2))

## ## A smooth plot of the age model can be achieved more simply via
## newdat3 <- with(dat, data.frame(Depth = seq(min(Depth), max(Depth),
##                                 length = 100)))
## pred3 <- transform(newdat3, Year = predict(mod, newdat3))
## ## plot...
## with(dat, plot(Depth, Date))
## with(pred3, lines(Depth, Year, col = "red"))
