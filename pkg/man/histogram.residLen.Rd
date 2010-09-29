\name{histogram.residLen}
\alias{histogram.residLen}
\title{Lattice histogram plot for residual lengths}
\description{
  Lattice \code{\link[lattice]{histogram}} method for
  \code{\link{residLen}} objects.
}
\usage{
\method{histogram}{residLen}(x, ..., xlab = NULL, ylab = NULL,
          type = c("percent", "count", "density"))
}
\arguments{
  \item{x}{Object of class \code{"residLen"}, the result of a call to
    \code{\link{residLen}}.}
  \item{xlab, ylab}{Axis labels. If not supplied, suitable defaults are
    generated, depending on whether RDA or CCA was used as the
    underlying ordination model.}
  \item{type}{Character string indicating type of histogram to be
    drawn. \code{"percent"} and \code{"count"} give relative frequency
    and frequency histograms, and can be misleading when breakpoints are
    not equally spaced. \code{"density"} produces a density scale
    histogram.

    See \code{\link[lattice]{histogram}} for further details.}
  \item{\dots}{Additional arguments passed to
    \code{\link[lattice]{histogram}}.}
}
\value{
  Returns an object of class \code{"trellis"}. See
  \code{\link[lattice]{histogram}} for details.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{residLen}}, \code{\link{plot.residLen}},
  \code{\link{hist.residLen}}, \code{\link{densityplot.residLen}}.}
\examples{
## load the Imbrie and Kipp example data
data(ImbrieKipp, SumSST, V12.122)

## squared residual lengths for Core V12.122
rlens <- residLen(ImbrieKipp, SumSST, V12.122)
rlens

## plot a histogram of the residual distances
histogram(rlens)

}
\keyword{hplot}
\keyword{methods}
