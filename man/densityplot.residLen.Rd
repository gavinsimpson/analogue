\name{densityplot.residLen}
\alias{densityplot.residLen}
\title{Lattice density plot for residual lengths}
\description{
  Lattice \code{\link[lattice]{densityplot}} method for
  \code{\link{residLen}} objects.
}
\usage{
\method{densityplot}{residLen}(x, ..., xlab = NULL, ylab = NULL)
}
\arguments{
  \item{x}{Object of class \code{"residLen"}, the result of a call to
    \code{\link{residLen}}.}
  \item{xlab, ylab}{Axis labels. If not supplied, suitable defaults are
    generated, depending on whether RDA or CCA was used as the
    underlying ordination model.}
  \item{\dots}{Additional arguments passed to
    \code{\link[lattice]{densityplot}}.}
}
\value{
  Returns an object of class \code{"trellis"}. See
  \code{\link[lattice]{densityplot}} for details.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{residLen}}, \code{\link{plot.residLen}},
  \code{\link{hist.residLen}}, \code{\link{histogram.residLen}}.}
\examples{
## load the Imbrie and Kipp example data
data(ImbrieKipp, SumSST, V12.122)

## squared residual lengths for Core V12.122
rlens <- residLen(ImbrieKipp, SumSST, V12.122)
rlens

## plot the density functions of the residual distances
densityplot(rlens)

}
\keyword{hplot}
\keyword{methods}
