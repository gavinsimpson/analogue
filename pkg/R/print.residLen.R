## print method
`print.residLen` <- function(x,
                             digits = min(3, getOption("digits") - 4),
                             probs = c(0.5, 0.75, 0.9, 0.95, 0.99), ...) {
    cat("\n")
    writeLines(strwrap("Squared residual lengths",
        prefix = "\t"))
    cat("\nCall:\n")
    cat(paste(deparse(x$call), "\n\n"))
    cat("Quantiles of training set lengths:\n")
    print(with(x, quantile(train, probs = probs)), digits = digits)
    cat("\nQuantiles of passive sample lengths:\n")
    print(with(x, quantile(passive, probs = probs)), digits = digits)
    cat("\n")
}
