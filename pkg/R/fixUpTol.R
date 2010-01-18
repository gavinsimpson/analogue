## fix-up tolerances for use in TF computations
fixUpTol <- function(tol, na.tol, small.tol, min.tol, f, env) {
    ## first missing tolerances
    if(any(NA.TOL <- is.na(tol) | is.infinite(tol))) {
        tol[NA.TOL] <-
            switch(na.tol,
                   min = min(tol, na.rm = TRUE),
                   mean = mean(tol, na.rm = TRUE),
                   max = max(tol, na.rm = TRUE))
    }
    ## second, replace tol < min.tol
    if(!is.null(min.tol) && any(MIN.TOL <- tol < min.tol)) {
        ## min.tol must be in or on extremesof range(env)
        if(min.tol < min(env) || min.tol > max(env))
            stop("'min.tol' must be >= min(env) and <= max(env)")
        if(small.tol == "fraction")
            frac <- f * diff(range(env))
        tol[MIN.TOL] <-
            switch(small.tol,
                   fraction = frac,
                   absolute = min.tol,
                   min = min(tol[tol >= min.tol], na.rm = TRUE))
    }
    return(tol)
}
