## Test `Stratiplot()` function

## load packages
library("testthat")
library("analogue")

context("Testing Stratiplot()")

test_that("Stratiplot() handles a vector of user-supplied labels", {
    ## Custom variable labels using expressions
    df <- data.frame(Age = 1:10, Var1 = rnorm(10), Var2 = rnorm(10),
                     Var3 = rnorm(10))
    ## Use a vector of expressions to label variables on plot
    ## See ?plotmath for syntax of expressions
    exprs <- expression(delta^{15}*N,       # label for Var1
                        delta^{18}*O,       # label for Var2
                        delta^{13}*C)       # label for Var3
    expect_silent(Stratiplot(Age ~ ., data = df, labelValues = exprs,
                             varTypes = "absolute"))
    ## Check the correct expr is in the environment of plt$axis()
    plt <- Stratiplot(Age ~ ., data = df, labelValues = exprs,
                      varTypes = "absolute")
    expect_identical(exprs, get("userLabels", envir = environment(plt$axis)))
})
