\name{hist.residLen}
\alias{hist.residLen}
\title{Histogram plot for residual lengths}
\description{
  Base graphics histogram plot method for \code{\link{residLen}}
  objects.
}
\usage{
\method{hist}{residLen}(x, breaks = "Sturges", freq = TRUE,
     probs = c(0.9, 0.95, 0.99), ncol = 1, lcol = "red",
     llty = "dashed", xlab = NULL, ylab = NULL,
     main = "Residual distances", rug = TRUE, ...)
}
\arguments{
  \item{x}{Object of class \code{"residLen"}, the result of a call to
    \code{\link{residLen}}.}
  \item{breaks}{How breakpoints for the histogram are determined. See
    \code{hist} for more details.}
  \item{freq}{logical; if \code{TRUE}, the histogram graphic is a
    representation of frequencies, the \code{counts} component of the
    result; if \code{FALSE}, probability densities, component
    \code{density}, are plotted (so that the histogram has a total area
    of one). Defaults to \code{TRUE}}
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
  \item{\dots}{additional arguments passed to \code{hist}.}
}
\value{
  A plot on the current device.
  
  Returns a list with two components (\code{train} and \code{passive}),
  each of which is an object returned by \code{hist}.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{residLen}}, \code{\link{plot.residLen}},
  \code{\link{histogram.residLen}}, \code{\link{densityplot.residLen}}.}
\examples{
data(swapdiat, swappH, rlgh)

## squared residual lengths for RLGH
rlens <- residLen(swapdiat, swappH, rlgh)
rlens

## plot a histogram of the residual distances
hist(rlens)

}
\keyword{hplot}
\keyword{methods}
