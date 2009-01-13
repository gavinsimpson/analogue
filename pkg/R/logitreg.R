`logitreg` <- function(object, groups, k = 1, ...)
    UseMethod("logitreg")

`logitreg.default` <- function(object, groups, k = 1, ...) {
    if(!is.factor(groups))
        groups <- factor(groups)
    lev <- levels(groups)
    #n.g <- length(lev)
    within <- without <- vector(mode = "list", length = length(lev))
    names(within) <- names(without) <- lev
    models <- vector(mode = "list", length = length(lev) + 1)
    names(models) <- c(lev, "Combined")
    k <- seq_len(k) + 1
    for(l in lev) {
        inds <- groups == l
        IN <- as.numeric(apply(object[inds, inds], 2,
                               function(x, k) {x[order(x)[k]]}, k = k))
        OUT <- as.numeric(apply(object[inds, !inds], 2,
                                function(x, k) {x[order(x)[k]]}, k = k))
        analogs <- rep(c(TRUE, FALSE), times = c(length(IN), length(OUT)))
        Dij <- c(IN, OUT)
        #dat <- data.frame(analogs, Dij)
        models[[l]] <- glm(analogs ~ Dij, data = data.frame(analogs, Dij),
                           family = binomial(link = "logit"))
        models[[l]]$Dij <- Dij
        within[[l]] <- IN
        without[[l]] <- OUT
    }
    IN <- do.call(c, within)
    OUT <- do.call(c, without)
    analogs <- rep(c(TRUE, FALSE), times = c(length(IN), length(OUT)))
    Dij <- c(IN, OUT)
    #dat <- data.frame(analogs, Dij)
    models[["Combined"]] <- glm(analogs ~ Dij,
                                data = data.frame(analogs, Dij),
                                family = binomial(link = "logit"))
    models[["Combined"]]$Dij <- Dij
    #models <- lapply(models, function(x) {class(x) <- "glm"; x})
    class(models) <- "logitreg"
    if(!is.null(attr(object, "method")))
        attr(models, "method") <- attr(object, "method")
    return(models)
}

print.logitreg <- function(x, ...) {
    nams <- names(x)
    N <- length(x)
    cat("\n")
    writeLines(strwrap("Object of class: \"logitreg\""))
    writeLines(strwrap(paste("Number of models:", N)))
    cat("\n")
    writeLines(strwrap("For groups:"))
    print(nams, ...)
    cat("\n")
    invisible(x)
}
