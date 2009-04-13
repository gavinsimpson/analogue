## compute the squared residual length statistic for a constrained
## ordination given a single constraint

`residLen` <- function(train, ...) {
    UseMethod("residLen")
}

`residLen.default` <- function(train, env, passive,
                               method = c("cca", "rda"), ...) {
    ## inline functions
    ## fitted values for passive samples
    fittedSuppl <- function(ord, newdata, csum) {
        #########################################################
        ## Fitted values for passive samples
        ## Arguments:
        ## ord     = vegan ordination object (CCA/RDA)
        ## newdata = matrix of passive species data with same
        ##           columns as that used to fit ord
        ## csum    = colum sums for the training data
        #########################################################
        ## species scores
        b <- predict(ord, type = "sp")
        ## site scores
        xi <- predict(ord, newdata = newdata, type = "wa")
        fik <- sweep(1 + (xi %*% t(b)), 2, csum / sum(csum),
                     "*")
        ## fitted values
        fik <- sweep(fik, 1, rowSums(newdata), "*")
    }
    sqrl.uni <- function(X, colsum, fik) {
        #########################################################
        ## Squared residual length for unimodal methods
        ## Arguments:
        ## X      = species data (training or passive) to which
        ##          we want to find the residual length
        ## colsum = column sums for the training species data
        ## fik    = fitted species data for either the training
        ##          or passive samples
        #########################################################
        yip <- rowSums(X)
        ypk <- colsum
        ypp <- sum(ypk)
        A <- (sweep(X, 2, ypk, "/") - sweep(fik, 2, ypk, "/"))^2
        B <- sweep(A, 2, ypk, "*")
        C <- rowSums(B / ypp)
        res <- (1 / (yip / ypp))^2 * C
        return(res)
    }
    ## merge train and passive
    dat <- join(train, passive)
    train <- dat[[1]]
    passive <- dat[[2]]
    ## check env is same length as nrow(train)
    if(!isTRUE(all.equal(length(env), nrow(train))))
        stop("'train' and 'env' imply different numbers of observations")
    ## ordinate
    if(missing(method))
        method <- "cca"
    method <- match.arg(method)
    ##if(method == "rda")
    ##    .NotYetUsed("method = \"rda\" not currently implemented.")
    FUN <- match.fun(method)
    ## ordinate
    ord <- FUN(X = train, Y = env)
    ## fitted values
    fit <- fitted(ord, type = "response")
    ## colSums of train
    ypk <- colSums(train)
    ## predict locations for the passive
    pred <- fittedSuppl(ord, passive, ypk)
    ## system.call
    .call <- match.call()
    .call[[1]] <- as.name("residLen")
    ## residual lengths for train
    res <- list(train = sqrl.uni(train, ypk, fit),
                passive = sqrl.uni(passive, ypk, pred),
                ordination = ord,
                call = .call)
    class(res) <- "residLen"
    attr(res, "method") <- method
    return(res)
}
