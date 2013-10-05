## Tests for the new distance compiled code

## load packages
library("testthat")
library_if_available("analogue")

context("Testing new distance compiled code functions")

## simple example using dummy data
train <- data.frame(matrix(abs(runif(200)), ncol = 10))
rownames(train) <- LETTERS[1:20]
colnames(train) <- as.character(1:10)
fossil <- data.frame(matrix(abs(runif(100)), ncol = 10))
colnames(fossil) <- as.character(1:10)
rownames(fossil) <- letters[1:10]

## test methods for x and y
test_that("distance matches compiled versions for x and y", {

    ## default settings
    expect_equal(distance(train, fossil),
                 oldDistance(train, fossil))

    ## euclidean
    method <- "euclidean"
    expect_equal(distance(train, fossil, method = method),
                 oldDistance(train, fossil, method = method))

})
