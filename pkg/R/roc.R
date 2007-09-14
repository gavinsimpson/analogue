roc <- function(object, groups, ...) UseMethod("roc")

roc.default <- function(object, groups, ...) {
  stop("Default roc method not yet implemented")
}

roc.mat <- function(object, groups, ...) {
  #######################################################################
  ## Calculation of the ROC curve points themselves,
  ## based on slightly modified code of Thomas Lumley
  ## (used with permission), in his Programmer's Niche article in R News
  ## (Vol. 4(1) 33--36). Uses the optimisations in the article to
  ## calculate the ROC curve itself.
  #######################################################################
  ## dists == T in Lumley's ROC
  #dists <- round(object$Dij[lower.tri(object$Dij)],2)
  dists <- object$Dij[lower.tri(object$Dij)]
  ## grps == D in Lumley's ROC
  grps <- outer(groups, groups, "==")
  grps <- grps[lower.tri(grps)]
  cutpoints <- rev(sort(unique(dists)))
  ## here, for our data, negating the first argument below (dist),
  ## produces a ROC curve that extends down, and to the right, below the
  ## 1:1 line of no discrimative power. So we do not use "-dist" below
  ## and the curve now moves up and to the left, as conventional. This
  ## does not change the ROC curve or it's interpretation.
  ## From Lumley(2004), would be "tab <- table(-dists, grps)"
  tab <- table(dists, grps)
  ## TPF == sens
  TPF <- cumsum(tab[,2]) / sum(tab[,2])
  ## FPE == mspec
  FPE <- cumsum(tab[,1]) / sum(tab[,1])
  roc.values <- TPF - FPE
  optimal <- as.numeric(names(which.max(roc.values)))
  names(FPE) <- names(TPF) <- names(roc.values) <- NULL
  wilcox <- wilcox.test(dists[grps == FALSE], dists[grps == TRUE],
                        conf.int = TRUE)
  n.within <- sum(grps)
  n.without <- length(dists) - n.within
  AUC <- wilcox$statistic / (n.without * n.within)
  retval <- list(TPF = TPF, FPE = FPE, roc.points = cutpoints,
                 roc.values = roc.values, optimal = optimal,
                 wilcox = wilcox, AUC = AUC,
                 n.within = n.within, n.without = n.without,
                 group = grps, dissims = dists,
                 method = object$method, call = sys.call(),
                 tab = tab)
  class(retval)<-"roc"
  retval
}

roc.analog <- function(object, groups, ...) {
  #######################################################################
  ## Calculation of the ROC curve points themselves,
  ## based on slightly modified code of Thomas Lumley
  ## (used with permission), in his Programmer's Niche article in R News
  ## (Vol. 4(1) 33--36). Uses the optimisations in the article to
  ## calculate the ROC curve itself.
  #######################################################################
  if(is.null(object$train))
    stop("'object$train' missing. Refit 'object' with argument 'keep.train = TRUE'")
  ## dists == T in Lumley's ROC
  #dists <- round(object$Dij[lower.tri(object$Dij)],2)
  dists <- object$train[lower.tri(object$train)]
  ## grps == D in Lumley's ROC
  grps <- outer(groups, groups, "==")
  grps <- grps[lower.tri(grps)]
  cutpoints <- rev(sort(unique(dists)))
  ## here, for our data, negating the first argument below (dist),
  ## produces a ROC curve that extends down, and to the right, below the
  ## 1:1 line of no discrimative power. So we do not use "-dist" below
  ## and the curve now moves up and to the left, as conventional. This
  ## does not change the ROC curve or it's interpretation.
  ## From Lumley(2004), would be "tab <- table(-dists, grps)"
  tab <- table(dists, grps)
  ## TPF == sens
  TPF <- cumsum(tab[,2]) / sum(tab[,2])
  ## FPE == mspec
  FPE <- cumsum(tab[,1]) / sum(tab[,1])
  roc.values <- TPF - FPE
  optimal <- as.numeric(names(which.max(roc.values)))
  names(FPE) <- names(TPF) <- names(roc.values) <- NULL
  wilcox <- wilcox.test(dists[grps == FALSE], dists[grps == TRUE],
                        conf.int = TRUE)
  n.within <- sum(grps)
  n.without <- length(dists) - n.within
  AUC <- wilcox$statistic / (n.without * n.within)
  .call <- match.call()
  .call[[1]] <- as.name("roc")
  retval <- list(TPF = TPF, FPE = FPE, roc.points = cutpoints,
                 roc.values = roc.values, optimal = optimal,
                 wilcox = wilcox, AUC = AUC,
                 n.within = n.within, n.without = n.without,
                 group = grps, dissims = dists,
                 method = object$method, call = .call,
                 tab = tab)
  class(retval)<-"roc"
  retval
}

print.roc <- function(x, digits = min(3, getOption("digits") - 4),
                         ...) {
  cat("\n")
  writeLines(strwrap("ROC curve of dissimilarities", prefix = "\t"))
  cat("\n")
  cat(paste("Optimal Dissimilarity =", round(x$optimal, digits), "\n\n"))
  cat(paste("AUC = ", round(x$AUC, digits),
            ", p-value: ", format.pval(x$wilcox$p.value), "\n", sep = ""))
  cat(paste("No. within:", x$n.within,
            "  No. outside:", x$n.without, "\n"))
  cat("\n")
  invisible(x)
}
