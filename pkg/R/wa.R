`wa` <- function(x, ...) UseMethod("wa")

`wa.default` <- function(x, env,
                         deshrink = c("inverse", "classical", "expanded", "none"),
                         tol.dw = FALSE, useN2 = TRUE,
                         na.tol = c("min","mean","max","wacalib"),
                         small.tol = c("min","fraction","absolute"),
                         min.tol = NULL, f = 0.1, ...)
{
    ## tol.dw is not yet implemented
    ##if (tol.dw)
    ##    .NotYetUsed("tol.dw", error = FALSE)
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
    ## first missing tolerances
    if(any(NA.TOL <- is.na(tol))) {
        tol[NA.TOL] <- switch(na.tol,
                              min = min(tol, na.rm = TRUE),
                              mean = mean(tol, na.rm = TRUE),
                              max = max(tol, na.rm = TRUE),
                              wacalib = .NotYetUsed("na.tol = \"wacalib\""))
    }
    ## second, replace tol < min.tol
    if(!is.null(min.tol) && any(MIN.TOL <- tol < min.tol)) {
        if(missing(small.tol))
            small.tol <- "min"
        small.tol <- match.arg(small.tol)
        ## min.tol must be in or on extremesof range(env)
        if(min.tol < min(env) || min.tol > max(env))
            stop("'min.tol' must be >= min(env) and <= max(env)")
        ## warn if "fraction" and f * diff(range(env))) < min.tol
        if(small.tol == "fraction") {
            if(!(f > 0 && f < 1))
                stop("'f' must be 0 < f < 1")
            frac <- f * diff(range(env))
            if(frac < min.tol)
                warning("Requested fraction of gradient is < minimum tolerance.")
        }
        tol[MIN.TOL] <- switch(small.tol,
                               fraction = frac,
                               absolute = min.tol,
                               min = min(tol[tol >= min.tol], na.rm = TRUE))
    }
    ## calculate WA estimate of env for each site
    if(tol.dw) {
        tol2 <- tol^2
        wa.env <- rowSums(sweep(sweep(x, 2, wa.optima, "*",
                                      check.margin = FALSE),
                                2, tol2, "/", check.margin = FALSE)) /
                                    rowSums(sweep(x, 2, tol2, "/",
                                                  check.margin = FALSE))
    } else {
        wa.env <- ((x %*% wa.optima) / rowSums(x))[,1, drop = TRUE]
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
                orig.x = x, orig.env = env)
    class(res) <- "wa"
    return(res)
}
