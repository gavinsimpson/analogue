## plot a principle curve in PCA space
plot.prcurve <- function(x, data, axes = 1:2,
                         seg = TRUE,
                         col.seg = "forestgreen",
                         col.curve = "red",
                         lwd.curve = 2, ...) {
    scl <- 0
    ordi <- rda(data)
    pred <- predict(ordi, x$s, type = "wa", scaling = scl)[,axes]
    scrs <- scores(ordi, display = "sites", scaling = scl,
                   choices = axes)
    xlim <- range(scrs[,1], pred[,1])
    ylim <- range(scrs[,2], pred[,2])
    plot(ordi, display = "sites", scaling = scl, type = "n",
         xlim = xlim, ylim = ylim, choices = axes, ...)
    points(scrs, ...)
    if(seg)
        segments(scrs[,1], scrs[,2], pred[,1], pred[,2],
                 col = col.seg)
    lines(pred[x$tag, 1:2], lwd = lwd.curve, col = col.curve)
}
