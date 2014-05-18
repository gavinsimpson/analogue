\name{varExpl}
\alias{varExpl}
\alias{varExpl.default}
\alias{varExpl.cca}
\alias{varExpl.prcurve}

\title{
  Variance explained by ordination axes
}
\description{
  Extracts information about the variance explained by ordination axes
  and expresses it in a variety of ways.
}
\usage{
varExpl(object, \dots)

\method{varExpl}{cca}(object, axes = 1L, cumulative = FALSE,
        pcent = FALSE, \dots)

\method{varExpl}{prcurve}(object, pcent = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object of an appropriate type. Currently only for
    objects that inherit from classes \code{"cca"} or \code{"prcurve"}.}
  \item{axes}{numeric vector indicating which axes to compute variance
    explained for.}
  \item{cumulative}{logical; should the variance be explained as a
    cumulative sum over the axes?}
  \item{pcent}{logical; should the variance explained be expressed as a
    percentage of the total variance.}
  \item{\dots}{additional arguments passed to other methods. Currently
    not used.}
}
%\details{
%
%}
\value{
  A numeric vector variance explained by each axis.
}
%\references{
%
%}
\author{
  Gavin L. Simpson
}

%\note{
%
%}

\seealso{
  See \code{\link{cca}} and \code{\link{prcurve}} for functions that
  produce objects that \code{varExpl()} can work with.
}
\examples{

data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit PCA
aber.pca <- rda(abernethy2)

## Distance along the first PCA axis
varExpl(aber.pca)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{utility}
