plot.roc <- function(x,
                     which = c(1:3,5),
                     prior = NULL,
                     show.stats = TRUE,
                     abline.col = "grey",
                     inGroup.col = "red",
                     outGroup.col = "blue",
                     caption = c("ROC curve",
                       "Dissimilarity profiles",
                       "TPF - FPF vs Dissimilarity",
                       "Likelihood ratios",
                       "Pr (A+ | d)"),
                     legend = "topright",
                     ask = prod(par("mfcol")) < length(which) &&
                     dev.interactive(),
                     ...) {
  if(!inherits(x, "roc"))
    stop("Plot method only for objects of class \"roc\".")
  if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
    stop("'which' must be in 1:5")
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if (ask) {
    par(ask = TRUE)
  }
  if(any(show[4:5])){
    ## removed argument which from bayesF to see if it is really needed
    #l.ratios <- bayesF(x, which = "both", prior = prior)
    l.ratios <- bayesF(x, prior = prior)
  }
  if(show[1]) {
    plot(x$FPE, x$TPF, type = "n",
         ylab = "TPF (sensitivity)",
         xlab = "1 - TNF (1 - specificity)")
    lines(x$FPE, x$TPF, ...)
    abline(0, 1, col = abline.col)
    mtext(caption[1], side = 3, line = 1.7, font = 2)#,
          #cex = par("cex.main"))
    if(show.stats) {
      txt <- paste("AUC =", round(x$AUC, 3))
      text(x = 0.9, y = 0, labels = txt, pos = 2, cex = 0.8) # was pos 3
    }
  }
  if(show[2]) {
    dens.in <- density(x$dissims[x$group == 1])
    dens.out <- density(x$dissims[x$group == 0])
    xlims <- switch(x$method,
                    SQchord = c(0,2),
                    chord = c(0, sqrt(2)),
                    c(0, max(x$dissims)))
    ylims <- range(0, dens.in$y, dens.out$y)
    plot(dens.in$x, dens.in$y, type = "n", axes = FALSE,
         xlim = xlims, ylim = ylims,
         ylab = "Density",
         xlab = paste("Dissimilarity (", x$method, ")", sep = ""))
    abline(h = 0, col = abline.col)
    lines(dens.in, col = inGroup.col)
    lines(dens.out, col = outGroup.col)
    abline(v = x$optimal, lty = "dotted", col = abline.col)
    axis(side = 2)
    axis(side = 1)
    box()
    mtext(caption[2], side = 3, line = 1.7, font = 2)#,
          #cex = par("cex.main"))
    legend("topright", legend = c("Analogue", "Not Analogue"),
           col = c(inGroup.col, outGroup.col),
           lty = "solid", lwd = 1,
           bty = "n", cex = 0.8,
           inset = 0.01)
  }
  if(show[3]) {
    cutpoints <- rev(x$roc.points)
    plot(cutpoints, x$roc.values, type = "n",
         ylab = "TPF - (1 - TNF)",
         xlab = paste("Dissimilarity (", x$method, ")", sep = ""))
    abline(h = 0, col = abline.col)
    lines(cutpoints, x$roc.values)
    abline(v = x$optimal, lty = "dotted", col = abline.col)
    mtext(caption[3], side = 3, line = 1.7, font = 2)#,
          #cex = par("cex.main"))
  }
  if(show[4]) {
    #margins <- par("mar")
    #par(mar = margins + c(0,0,0,2))
    #l.ratios <- bayesF(x, which = "both")
    #spec <- 1 - x$FPE
    dissims <- rev(x$roc.points)
    plot(dissims, l.ratios$pos, type = "n", axes = FALSE,
         ylab = "LR (+)",
         xlab = paste("Dissimilarity (", x$method, ")", sep = ""))
    abline(v = x$optimal, lty = "dotted", col = abline.col)
    lines(dissims, l.ratios$pos, col = inGroup.col)
    axis(side = 1)
    axis(side = 2)
    usr <- par("usr")
    finite.vals <- is.finite(l.ratios$neg)
    rany <- (max(l.ratios$neg[finite.vals]) -
             min(l.ratios$neg[finite.vals])) * 0.04
    par(usr = c(usr[1:2], min(l.ratios$neg[finite.vals]) - rany,
          max(l.ratios$neg[finite.vals]) + rany))
    lines(dissims, l.ratios$neg, col = outGroup.col)
    axis(side = 4)
    #mtext("LR (-)", side = 4, line = 3)
    box()
    mtext(caption[4], side = 3, line = 1.7, font = 2)
    legend("topright", legend = c("LR (+)", "LR (-)"),
           col = c(inGroup.col, outGroup.col),
           lty = "solid", lwd = 1,
           bty = "n", cex = 0.8,
           inset = 0.01)
  }
  if(show[5]) {
    #prob.pos <- l.ratios$posterior$pos / (1 + x$posterior$pos)
    #plot(rev(x$roc.points), prob.pos, type = "l", col = "red")
    plot(l.ratios, abline.col = abline.col, col = inGroup.col,
         ylab = "Pr (A+ | d)",
         xlab = paste("Dissimilarity (", x$method, ")", sep = ""))
    mtext(caption[5], side = 3, line = 1.7, font = 2)
  }
  invisible()
}
