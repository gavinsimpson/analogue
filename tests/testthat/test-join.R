## Test `join()` function

## load packages
library("testthat")
library("analogue")

context("Testing join()")

## Use example from #16 as provided by Richard Telford
spp <- data.frame(cat = 1:2, dog = 2:3, fish = 3:4)
fos <- data.frame(mouse = 7:8, dog = 9:10, cat = 5:6)
row.names(spp) <- row.names(fos) <- as.character(1:2)
LEFT <- structure(list(spp = spp,
                       fos = cbind(fos[, c("cat", "dog")], fish = c(0,0))),
                  class = "join", type = "left")
OUTER <- structure(list(spp = cbind(spp, mouse = c(0,0)),
                        fos = cbind(fos[, c("cat", "dog")], fish = c(0,0),
                        mouse = fos[, "mouse"])),
                   class = "join", type = "outer")
cn <- intersect(colnames(spp), colnames(fos))
INNER <- structure(list(spp = spp[, cn], fos = fos[, cn]),
                   class = "join", type = "inner")

test_that("join() with type = 'left' works", {
    obj <- join(spp, fos, type = "left")
    expect_equal(LEFT, obj)
})

test_that("join() with type = 'outer' works", {
    obj <- join(spp, fos, type = "outer")
    expect_equal(OUTER, obj)
})

test_that("join() with type = 'inner' works", {
    obj <- join(spp, fos, type = "inner")
    expect_equal(INNER, obj)
})
