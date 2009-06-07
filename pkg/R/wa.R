`wa` <- function(x, ...) UseMethod("wa")

`wa.default` <- function(x, env,
                         deshrink = c("inverse", "classical", "expanded", "none"),
                         tol.dw = FALSE, useN2 = TRUE,
                         na.tol = c("min","mean","max"),
                         small.tol = c("min","fraction","absolute"),
                         min.tol = NULL, f = 0.1, ...)
{
    ## x = species abundances (weights), env = response vector
    x <- as.matrix(x)
    env <- as.numeric(env)
    ## drop species with no information
    if(any(csum <- colSums(x) == 0))
        x <- x[, !csum, drop = FALSE]
    ## sample summaries
    n.samp <- nrow(x)
    n.spp <- ncol(x)
    if(missing(deshrink))
        deshrink <- "inverse"
    deshrink <- match.arg(deshrink)
    ## calculate WA optima for each species in x
    wa.optima <- w.avg(x, env)
    ## compute tolerances
    tolerances <- tol <- w.tol(x, env, wa.optima, useN2 = useN2)
    ## fix-up tolerances for use in TF computations
    if(missing(na.tol))
        na.tol <- "min"
    na.tol <- match.arg(na.tol)
    if(missing(small.tol))
        small.tol <- "min"
    small.tol <- match.arg(small.tol)
    if(small.tol == "fraction") {
        if(!(f > 0 && f < 1))
            stop("'f' must be 0 < f < 1")
        frac <- f * diff(range(env))
        if(frac < min.tol)
            warning("Requested fraction of gradient is < minimum tolerance.")
    }
    tol <- fixUpTol(tol, na.tol = na.tol, small.tol = small.tol,
                    min.tol = min.tol, f = f, env = env)
    ## calculate WA estimate of env for each site
    if(tol.dw) {
        wa.env <- WATpred(x, wa.optima, tol)
    } else {
        wa.env <- WApred(x, wa.optima)
    }
    ## taken averages twice so deshrink
    expanded <- switch(deshrink,
                       inverse = inv.deshrink(env, wa.env),
                       classical = class.deshrink(env, wa.env),
                       expanded = expand.deshrink(env, wa.env),
                       none = no.deshrink(env, wa.env))
    wa.env <- expanded$env
    coefficients <- coef(expanded)
    ## site/sample names need to be reapplied
    names(wa.env) <- rownames(x)
    ## species names need to be reapplied
    names(wa.optima) <- colnames(x)
    ## residuals
    resi <- wa.env - env
    ## RMSE of predicted/fitted values
    rmse <- sqrt(mean(resi^2))
    ## r-squared
    r.squared <- cor(wa.env, env)^2
    ## bias statistics
    avg.bias <- mean(resi)
    max.bias <- maxBias(resi, env)
    ## the function call
    .call <- match.call()
    ## need to reset due to method dispatch
    .call[[1]] <- as.name("wa")
    ## returned object
    res <- list(wa.optima = wa.optima,
                tolerances = tolerances,
                model.tol = tol,
                fitted.values = wa.env,
                residuals = resi,
                coefficients = coefficients,
                rmse = rmse, r.squared = r.squared,
                avg.bias = avg.bias, max.bias = max.bias,
                n.samp = n.samp, n.spp = n.spp,
                deshrink = deshrink, tol.dw = tol.dw,
                call = .call,
                orig.x = x, orig.env = env,
                options.tol =
                list(useN2 = useN2,
                     na.tol = na.tol,
                     small.tol = small.tol,
                     min.tol = min.tol,
                     f = f))
    class(res) <- "wa"
    return(res)
}
