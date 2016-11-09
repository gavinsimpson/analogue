## Test `join()` function

## load packages
library("testthat")
library("analogue")

context("Testing join()")

## Use example from #16 as provided by Richard Telford
spp <- data.frame(cat = 1:2, dog = 2:3, fish = 3:4)
fos <- data.frame(mouse = 7:8, dog = 9:10, cat = 5:6)
row.names(spp) <- row.names(fos) <- as.character(1:2)
left <- structure(list(spp = spp,
                       fos = cbind(fos[, c("cat", "dog")], fish = c(0,0))),
                  class = "join", type = "left")

test_that("join() with type = 'left' works", {
    obj <- join(spp, fos, type = "left")
    expect_equal(left, obj)
})
