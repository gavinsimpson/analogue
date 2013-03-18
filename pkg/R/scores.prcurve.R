## `scores` extractor function for prcurve class
`scores.prcurve` <- function(x, ...) {
  matrix(x$lambda, ncol = 1)
}
