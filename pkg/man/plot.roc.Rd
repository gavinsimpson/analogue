\name{plot.roc}
\alias{plot.roc}
\concept{ROC}
\concept{likelihood ratios}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot ROC curves and associated diagnostics}
\description{
  Produces up to four plots (selectable by \code{"which"}) from the
  results of a call to \code{\link{roc}}, including the ROC curve
  itself.
}
\usage{
\method{plot}{roc}(x,
     which = c(1:3,5),
     group = "Combined",
     prior = NULL,
     show.stats = TRUE,
     abline.col = "grey",
     inGroup.col = "red",
     outGroup.col = "blue",
     caption = c("ROC curve", "Dissimilarity profiles",
                 "TPF - FPF vs Dissimilarity",
                 "Likelihood ratios"),
     legend = "topright",
     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"roc"}.}
  \item{which}{numeric vector; which aspects of \code{"roc"} object to
    plot if a subset of the plots is required, specify a subset of the
    numbers \code{1:5}.}
  \item{group}{character vector of length 1 giving the name of the group
        to plot.}
  \item{prior}{numeric vector of length 2 (e.g. \code{c(0.5, 0.5)})
    specifiying the prior probabilities of analogue and
    no-analogue. Used to generate posterior probability of analogue
    using Bayes factors in plot 5 (\code{which = 5}).}
  \item{show.stats}{logical; should concise summary statistics of the
    ROC analysis be displayed on the plot?}
  \item{abline.col}{character string or numeric value; the
    colour used to draw vertical lines at the optimal dissimilarity on
    the plots.}
  \item{inGroup.col}{character string or numeric value; the colour used
    to draw the density curve for the dissimilarities between sites in
    the same group.}
  \item{outGroup.col}{character string or numeric value; the colour used
    to draw the density curve for the dissimilarities between sites not in
    the same group.}
  \item{caption}{vector of character strings, containing the captions to
    appear above each plot.}
  \item{legend}{character; position of legends drawn on plots. See
    Details section in \code{\link{legend}} for keywords that can be
    specified.}
  \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, see \code{par(ask=.)}.}
  \item{\dots}{graphical arguments passed to other graphics functions.}
}
\details{
  This plotting function is modelled closely on \code{\link{plot.lm}}
  and many of the conventions and defaults for that function are
  replicated here.

  First, some definitions:
  
  \describe{
    \item{\strong{TPF}}{ True Positive Fraction, also known as
      \emph{sensitivity}.}
    \item{\strong{TNF}}{ True Negative Fraction, also known as
      \emph{specificity}.}
    \item{\strong{FPF}}{ False Positive Fraction; the complement of TNF,
      calculated as 1 - TNF. This is often referred to a 1 -
      specificity. A false positive is also known as a type I error.} 
    \item{\strong{FNF}}{ False Negative Fraction; the complement of TPF,
      calculated as 1 - TPF. A false negative is also known as a type II
      error. }
    \item{\strong{AUC}}{ The Area Under the ROC Curve.}
  }
  
  The "ROC curve" plot (\code{which = 1},) draws the ROC curve itself as
  a plot of the False Positive Fraction against the True Positive
  Fraction. A diagonal 1:1 line represents no ability for the
  dissimilarity coefficient to differentiate between groups. The AUC
  statistic may also be displayed (see argument \code{"show.stats"}
  above).

  The "Dissimilarity profile" plot (\code{which = 2}), draws the density
  functions of the dissimilarity values (\emph{d}) for the correctly
  assigned samples and the incorrectly assigned samples. A dissimilarity
  coefficient that is able to well distinguish the sample groupings will
  have density functions for the correctly and incorrectly assigned
  samples that have little overlap. Conversely, a poorly discriminating
  dissimilarity coefficient will have density profiles for the two
  assignments that overlap considerably. The point where the two curves
  cross is the optimal dissimilarity or critical value, \emph{d'}. This
  represents the point where the difference between TPF and FPF is
  maximal. The value of \emph{d} at the point where the difference
  between TPF and FPF is maximal will not neccesarily be the
  same as the value of \emph{d'} where the density profiles cross. This
  is because the ROC curve has been estimated at discrete points
  \emph{d}, which may not include excatly the optimal \emph{d'}, but
  which should be close to this value if the ROC curve is not sampled on
  too coarse an interval.

  The "TPF - FPF vs Dissimilarity" plot, draws the difference between
  the ROC curve and the 1:1 line. The point where the ROC curve is
  farthest from the 1:1 line is the point at which the ROC curve has
  maximal slope. This is the optimal value for \emph{d}, as discussed
  above.

  The "Likelihood ratios" plot, draws two definitions of the
  slope of the ROC curve as the likelihood functions LR(+), and
  LR(-). LR(+), is the likelihood ratio of a positive test result, that
  the value of \emph{d} assigns the sample to the group it belongs
  to. LR(-) is the likelihood ratio of a negative test result, that the
  value of \emph{d} assigns the sample to the wrong group.

  LR(+) is defined as \eqn{LR(+) = TPF / FPF} (or sensitivity / (1 -
  specificity)), and LR(-) is defined as \eqn{LR(-) = FPF / TNF} (or (1
  - sensitivity) / specificity), in Henderson (1993).

  The \dQuote{probability of analogue} plot, draws the posterior
  probability of analogue given a dissimilarity. This is the
  LR(+) likelihood ratio values multiplied by the prior odds of
  analogue, for given values of the dissimilarity, and is then converted
  to a  probability. 
}
\value{
  One or more plots, drawn on the current device.
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
\author{Gavin L. Simpson. Code borrows heavily from \code{\link{plot.lm}}.}
\seealso{\code{\link{roc}}}
\examples{
## continue the example from roc()
example(roc)

## draw the ROC curve
plot(swap.roc, 1)

## draw the four default diagnostic plots
opar <- par(mfrow = c(2,2))
plot(swap.roc)
par(opar)
}
\keyword{hplot}
\keyword{methods}
