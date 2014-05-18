\name{plot.residLen}
\alias{plot.residLen}
\title{Plot method for residual lengths}
\description{
  Base graphics plot method for \code{\link{residLen}} objects.
}
\usage{
\method{plot}{residLen}(x, probs = c(0.9, 0.95, 0.99), ncol = 1,
     lcol = "red", llty = "dashed", xlab = NULL, ylab = NULL,
     main = "Residual distances", rug = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"residLen"}, the result of a call to
    \code{\link{residLen}}.}
  \item{probs}{numeric; vector of probability quantiles to compute from
    the sets of residual distances.}
  \item{ncol}{numeric; number of columns for the plot layout. Choices
    are \code{1} or \code{2}. Determines whether the histograms are
    plotted above or beside each other.}
  \item{lcol, llty}{colour and line-type for the quantiles.}
  \item{xlab, ylab}{Axis labels. If not supplied, suitable defaults are
    generated, depending on whether RDA or CCA was used as the
    underlying ordination model.}
  \item{main}{character; title for the plot.}
  \item{rug}{logical; should rug plots of the actual distances be drawn?}
  \item{\dots}{additional arguments passed to \code{plot}.}
}
\value{
  A plot on the current device.
  
  Returns, invisibly, a list with two components (\code{train} and
  \code{passive}), each and object of the type returned by
  \code{density}.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{residLen}}, \code{\link{plot.residLen}},
  \code{\link{histogram.residLen}}, \code{\link{densityplot.residLen}}.}
\examples{
## load the Imbrie and Kipp example data
data(ImbrieKipp, SumSST, V12.122)

## squared residual lengths for Core V12.122
rlens <- residLen(ImbrieKipp, SumSST, V12.122)
rlens

## plot a histogram of the residual distances
plot(rlens)

}
\keyword{hplot}
\keyword{methods}
