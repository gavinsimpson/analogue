## compute Bayes factors (syn. likelihood ratios) of
## positive and negative events
bayesF <- function(x, prior = NULL) {
  if(class(x) != "roc")
    stop("Likelihood ratios only for objects of class 'roc'")
  ## only allow a single option to which
  #if(length(which) > 1)
  #  stop("'which' must be of a character string of length 1")
  ## which can be one of "both","pos","neg", check which is it
  #which <- pmatch(which, c("both","pos","neg"), nomatch = FALSE)
  #if(!which)
  #  stop(paste("Unknown option '", which, "' given for 'which'.",
  #             "\nMust be one of 'both', 'pos', or 'neg'.",
  #             sep = ""))
  n.obs <- length(x$FPE)
  #if(which == 1) {
    pos <- x$TPF / x$FPE
    neg <- (1 - x$TPF) / (1 - x$FPE)
  #} else if(which == 2) {
  #  pos <- x$TPF / x$FPE
  #  neg <- NULL
  #} else {
  #  neg <- (1 - x$TPF) / (1 - x$FPE)
  #  pos <- NULL
  #}
  tot.comp <- x$n.within + x$n.without
  ## prior probs
  if(is.null(prior)) {
    prior.prob.pos <- x$n.within / tot.comp
    prior.prob.neg <- x$n.without / tot.comp
  } else {
    prior.prob.pos <- prior[1]
    prior.prob.neg <- prior[2]
  }
  ## prior odds
  prior.odds.pos <- prior.prob.pos / (1 - prior.prob.pos)
  prior.odds.neg <- prior.prob.neg / (1 - prior.prob.neg)
  ## posterior odds
  post.odds.pos <- pos * prior.odds.pos
  post.odds.neg <- neg * prior.odds.neg
  ## return object
  retval <- list(pos = pos,
                 neg = neg,
                 posterior = list(pos = post.odds.pos,
                   neg = post.odds.neg),
                 prior = list(pos = prior.prob.pos,
                   neg = prior.prob.neg),
                 roc.points = x$roc.points,
                 optimal = x$optimal,
                 #method = x$method,
                 object = deparse(substitute(x))
                 )
  class(retval) <- "bayesF"
  return(retval)
}

print.bayesF <- function(x, digits = min(3, getOption("digits") - 4),
                         ...) {
  cat("\n")
  writeLines(strwrap("Bayes factors (likelihood ratios)", prefix = "\t"))
  cat("\n")
  cat(paste("Object:", x$object, "\n"))
  cat("\nPrior probabilities:\n")
  cat(paste("Positive:", round(x$prior$pos, digits),
            "  Negative:", round(x$prior$neg,digits), "\n"))
  cat("\n")
  invisible(x)
}
