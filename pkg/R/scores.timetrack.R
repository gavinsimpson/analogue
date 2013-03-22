`scores.timetrack` <- function(x, which = c("passive","ordination"),
                               scaling = x$scaling, choices = 1:2,
                               ...) {
  which <- match.arg(which)
  scrs <- if(which == "passive") {
    fitted(x, which = which, choices = choices, ...)
  } else {
    scores(x, ..., choices = choices, scaling = scaling)
  }
  scrs
}
