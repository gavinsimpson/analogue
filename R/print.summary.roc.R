`print.summary.roc` <- function(x, digits = 3,
                              print.gap = 2, ...) {
    cat("\n")
    writeLines(strwrap("ROC curves of dissimilarities:"))
    cat("\n")
    class(x) <- "data.frame"
    x$`p-value` <- format.pval(x$`p-value`)
    print(x, digits = digits, print.gap = print.gap, ...)
    invisible(x)
}
