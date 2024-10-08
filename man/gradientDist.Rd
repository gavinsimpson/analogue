\name{gradientDist}
\alias{gradientDist}
\alias{gradientDist.default}
\alias{gradientDist.cca}
\alias{gradientDist.prcurve}

\title{
  Positions of samples along a unit-length ordination gradient.
}
\description{
  Extracts information as to the locations of samples along an
  ordination gradient. \code{gradientDist()} standardises the entire
  gradient to the interval 0, \ldots, 1, to allow comparison between
  methods or data sets.
}
\usage{
gradientDist(object, \dots)

\method{gradientDist}{default}(object, na.rm = TRUE, \dots)

\method{gradientDist}{cca}(object, na.rm = TRUE, axis = 1L,
             scaling = 0, \dots)

\method{gradientDist}{prcurve}(object, na.rm = TRUE, \dots)
}

\arguments{
  \item{object}{an R object of an appropriate type. For the default
    method, any R object that can be coerced to a vector.}
  \item{na.rm}{logical; should missing values be removed?}
  \item{axis}{numeric, length 1; the ordination axis to take as the
    gradient.}
  \item{scaling}{Scaling to apply to the site scores. Default is to do
    no scaling. See \code{\link[vegan]{scores.cca}} for details.}
  \item{\dots}{additional arguments passed to other methods. In the
    \code{"cca"} method, these are also passed to
    \code{\link[vegan]{scores.cca}}.}
}
%\details{
%
%}
\value{
  A numeric vector of positions along the gradient, scaled to the range
  0, \ldots, 1.
}
%\references{
%
%}
\author{
  Gavin L. Simpson
}

\seealso{
  See \code{\link[vegan]{cca}} and \code{\link{prcurve}} for functions that
  produce objects that \code{gradientDist()} can work with.
}
\examples{

data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit PCA
aber.pca <- rda(abernethy2)

## Distance along the first PCA axis
gradientDist(aber.pca)
}

\keyword{multivariate}
\keyword{utility}
