lines.prcurve <- function(x, data, axes = 1:2, segments = TRUE,
                          col = "red", col.seg = "forestgreen",
                          lwd = 2, lwd.seg = 1,
                          ...) {
  scl <- 0
  ordi <- rda(data)
  pred <- predict(ordi, x$s, type = "wa", scaling = scl)[, axes]
  scrs <- scores(ordi, display = "sites", scaling = scl, choices = axes)
  if(segments)
    segments(scrs[, 1], scrs[, 2], pred[, 1], pred[, 2], 
             col = col.seg, lwd = lwd.seg)
  lines(pred[x$tag, 1:2], lwd = lwd, col = col,
        ...)
  invisible()
}
