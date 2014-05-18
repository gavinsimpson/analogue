\name{bayesF}
\alias{bayesF}
\alias{print.bayesF}
\alias{plot.bayesF}
\concept{likelihood ratios}
\concept{Bayes factors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayes factors}
\description{
  Calculates Bayes factors or likelihood ratios of analogue and
  no-analogue results.
}
\usage{
bayesF(x, prior = rep(0.5, 2))

\method{plot}{bayesF}(x, group = "all", xlab = NULL, ylab = "Pr (A+ | d)",
        col = "red", abline.col = "lightgrey", abline.lty = "dashed", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{for \code{bayesF} an object of class \code{roc}. For the plot
        method, an object of class \code{bayesF}, usually the result of a
        call to \code{bayesF}.}
  \item{prior}{numeric; the prior probabilities of analogue and
    no-analogue, provided as a vector of length 2 whose elements sum to
    1. If not provided, the function will use the relative occurences of
    analogue and no analogue situations used to evaluate the ROC curve.}
  \item{group}{character vector of length 1 giving the name of the group
        to plot, or \code{"all"} to plot all groups in \code{x}.}
  \item{xlab,ylab}{the x- and y-axis labels for the plot.}
  \item{col}{colour of the line used to draw the posterior probability.}
  \item{abline.col}{colour of the vertical line drawn to indicate the
    optimal dissimilarity determined from the ROC curve.}
  \item{abline.lty}{Line type for indicator of optimal ROC dissimilarity
    threshold. See \code{\link{par}} for the allowed line types.}
  \item{\dots}{other plot arguments passed to plotting
    functions. Currently ignored.}
}
\details{
  LR(+), is the likelihood ratio of a positive test result, that
  the value of \emph{d} assigns the sample to the group it belongs
  to. LR(-) is the likelihood ratio of a negative test result, that the
  value of \emph{d} assigns the sample to the wrong group.

  LR(+) is defined as \eqn{LR(+) = TPF / FPF} (or sensitivity / (1 -
  specificity)), and LR(-) is defined as \eqn{LR(-) = FPF / TNF} (or (1
  - sensitivity) / specificity), in Henderson (1993).

  The posterior probability of analogue given a dissimilarity is the 
  LR(+) likelihood ratio values multiplied by the prior odds of
  analogue, for given values of the dissimilarity, and is then converted
  to a  probability.
  
  The plotting function currently only draws the posterior probability
  of analogue based on the Bayes factor or likelihood ratio of a
  positive event (analogue).
}
\value{
  For \code{plot.bayesF} a plot on the currently active device.

  For \code{bayesF}, a list containing the results of computing Bayes
  factors for each group in \code{x}. Each component of this list is
  itself a list with the following components:

  \item{bayesF, posterior.odds, posterior.probs, prior.prob}{Bayes
    factors, posterior odds and probabilities and prior probabilities of
    true analogue and true non-analogue events. Each components is a list
    with two components; \code{pos} (for true analogue events) and
    \code{neg} (for true non-analogue events). The components of
    \code{prior.prob} are vectors of length 1, whilst components of the
    other lists are numeric vectors.}
  \item{roc.points}{numeric; the points at which the ROC curve was
    evaluated.}
  \item{optimal}{numeric; the optimal dissimilarity as assessed by the
    ROC curve.}
  \item{max.roc}{numeric; the position along the ROC curve at which the
    slope of the ROC curve is maximal. This is the index of this point
    on the curve, and can be used to extract the element of
    \code{bayesF}, \code{posterior.odds} and \code{posterior.probs} for
    the optimal dissimilarity.}
  %\item{object}{name of the object passed as argument \code{x}.}
}
\references{
  
  Brown, C.D., and Davis, H.T. (2006) Receiver operating characteristics
  curves and related decision measures: A tutorial. \emph{Chemometrics
    and Intelligent Laboratory Systems} \bold{80}, 24--38.
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.

  Henderson, A.R. (1993) Assessing test accuracy and its clinical
  consequences: a primer for receiver operating characteristic curve
  analysis. \emph{Annals of Clinical Biochemistry} \strong{30},
  834--846.

}
\author{Gavin L. Simpson}
\seealso{\code{\link{roc}} and \code{\link{plot.bayesF}}.}
\examples{
## load the example data
data(swapdiat, swappH, rlgh)

## merge training and test set on columns
dat <- join(swapdiat, rlgh, verbose = TRUE)

## extract the merged data sets and convert to proportions
swapdiat <- dat[[1]] / 100
rlgh <- dat[[2]] / 100

## fit an analogue matching (AM) model using the squared chord distance
## measure - need to keep the training set dissimilarities
swap.ana <- analog(swapdiat, rlgh, method = "SQchord",
                   keep.train = TRUE)

## fit the ROC curve to the SWAP diatom data using the AM results
## Generate a grouping for the SWAP lakes
clust <- hclust(as.dist(swap.ana$train), method = "ward")
grps <- cutree(clust, 12)

## fit the ROC curve
swap.roc <- roc(swap.ana, groups = grps)
swap.roc

## calculate the Bayes factors of analogue and no-analogue
## (uses observed probabilities of analogue/no-analogue
swap.bayes <- bayesF(swap.roc)
swap.bayes

## plot the probability of analogue
plot(swap.bayes)

\dontrun{
## calculate the Bayes factors of analogue and no-analogue
## with prior probabilities c(0.5, 0.05)
swap.bayes2 <- bayesF(swap.roc, prior = c(0.5, 0.05))
swap.bayes

## plot the probability of analogue
plot(swap.bayes2)
}
}
\keyword{univar}
\keyword{methods}
