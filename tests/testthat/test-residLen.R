library("testthat")
library("analogue")

context("Squared residual lengths of passive samples")

data(ImbrieKipp, SumSST, V12.122)

test_that("training rows passed as passive samples reach the axis minima", {
    X <- ImbrieKipp
    passive <- X[1:5, ]
    cs <- colSums(X)
    cm <- colMeans(X)
    p <- cs / sum(cs)

    for (method in c("cca", "rda")) {
        result <- residLen(X, SumSST, passive, method = method)
        b <- predict(result$ordination, type = "sp")[, 1]
        minimum <- vapply(seq_len(nrow(passive)), function(i) {
            y <- unlist(passive[i, ])
            optimize(function(x) {
                if (method == "cca") {
                    sum((y / sum(y) - (1 + x * b) * p)^2 / p)
                } else {
                    mean((y - (cm + x * b))^2)
                }
            }, if (method == "cca") c(-10, 10) else c(-1000, 1000))$objective
        }, numeric(1))
        expect_equal(unname(result$passive), minimum, tolerance = 1e-8)
        expect_true(all(result$passive <= result$train[1:5] + 1e-8))

        ## Training residuals still use the constrained fitted positions.
        fit <- fitted(result$ordination, type = "response")
        expected <- if (method == "cca") {
            rowSums(sweep((as.matrix(X) / rowSums(X) -
                          fit / rowSums(X))^2, 2, p, "/"))
        } else {
            rowMeans((as.matrix(X) - fit)^2)
        }
        expect_equal(result$train, expected, tolerance = 1e-8)
    }
})

test_that("passive fits are projections onto the constrained species subspace", {
    X <- ImbrieKipp
    env <- data.frame(SumSST, curved = SumSST^2)
    passive <- join(X, V12.122, type = "left")[[2]][1:4, ]
    ## Unequal totals exercise the return from CCA profiles to abundances.
    passive <- passive * c(0.5, 2, 3, 1)

    for (method in c("cca", "rda")) {
        for (naxes in 1:2) {
            result <- residLen(X, env[, seq_len(naxes), drop = FALSE],
                               passive, method = method)
            fit <- fitted(result$ordination, type = "response")
            for (npassive in c(1L, 4L)) {
                y <- as.matrix(passive[seq_len(npassive), , drop = FALSE])
                if (method == "cca") {
                    centre <- colSums(X) / sum(as.matrix(X))
                    weight <- sqrt(centre)
                    training <- sweep(fit / rowSums(X), 2, centre, "-")
                    response <- sweep(y / rowSums(y), 2, centre, "-")
                } else {
                    centre <- colMeans(X)
                    weight <- rep(1, ncol(X))
                    training <- sweep(fit, 2, centre, "-")
                    response <- sweep(y, 2, centre, "-")
                }
                training <- sweep(training, 2, weight, "/")
                response <- sweep(response, 2, weight, "/")
                ## Derive an orthonormal basis from fitted response data,
                ## independently of vegan's site scores and eigenvalues.
                basis <- svd(training, nu = 0, nv = naxes)$v
                projection <- response %*% basis %*% t(basis)
                expected <- rowSums((response - projection)^2)
                if (method == "rda") {
                    expected <- expected / ncol(X)
                }
                expect_equal(unname(result$passive[seq_len(npassive)]), unname(expected),
                             tolerance = 1e-8)
                reconstructed <- sweep(projection, 2, weight, "*")
                reconstructed <- sweep(reconstructed, 2, centre, "+")
                if (method == "cca") {
                    reconstructed <- reconstructed * rowSums(y)
                }
                dimnames(reconstructed) <- dimnames(y)
                expect_equal(analogue:::fittedY(result$ordination, y,
                                                colSums(X)),
                             reconstructed, tolerance = 1e-8)
            }
        }
    }
})

test_that("residLen() accepts a single passive sample", {
    for (method in c("cca", "rda")) {
        batch <- residLen(ImbrieKipp, SumSST, V12.122[1:2, ], method = method)
        single <- residLen(ImbrieKipp, SumSST, V12.122[1, , drop = FALSE],
                           method = method)
        expect_equal(single$passive, batch$passive[1])
        expect_equal(single$train, batch$train)
    }
})
