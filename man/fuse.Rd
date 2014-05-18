\name{fuse}
\alias{fuse}
\alias{fuse.matrix}
\alias{fuse.dist}
\title{Fused dissimilarities}
\description{
  Combines dissimilarities from two or more dissimilarity objects into a
  single dissimilarity object so that both original dissimilarities
  contribute equally. Weighted combinations of the original objects can
  also be created.
}
\usage{
fuse(..., weights = NULL)

\method{fuse}{matrix}(..., weights = NULL)

\method{fuse}{dist}(..., weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{objects to fuse. Methods currently exist for objects of
    class \code{"matrix"} and \code{"dist"} objects. Method dispatch is
    performed on the first object specified in \code{\dots}. A minimum
    of two objects must be supplied.}
  \item{weights}{numeric; vector of weights that sum to 1, one weight
    for each object supplied in \code{\dots}.}
}
\details{
  Fuses, or combines, dissimilarity objects in a very flexible way to
  create a single dissimilarity object that incorporates the separate
  dissimilarities. In analogue matching, we may wish to combine
  information from two or more proxies, such as diatoms and cladocera,
  or from biological and chemical or physical data in the case of
  matching modern samples.

  The function can also be used to fuse dissimilarity objects created
  from a single data set but using different dissimilarity
  coefficients. In this way one could create a new dissimilarity object
  combining dissimilarity based on abundance data and presence absence
  data into a single measure.

  \code{fuse} uses the method of Melssen et al. (2006) to combine
  dissimilarities. The dissimilarities in each dissimilarity object are
  scaled so that the maximum dissimilarity in each object is 1. The
  scaled dissimilarity objects are then weighted according to the
  supplied weights. If no weights are supplied (the default) the
  dissimilarity objects are weighted equally; \code{weights = rep(1/N,
    N)}, where \code{N} is the number of dissimilarity objects fused.
  
  \deqn{D_{fused}(j, k) = \sum_{i = 1}^N w_i D_{ijk}}{%
    D[fused](j, k) = sum(w[i] * D[ijk])}
  where \eqn{D_{fused}(j, k)}{D[fused](j, k)} is the fused dissimilarity
  between samples \eqn{j} and \eqn{k}, \eqn{w_i}{w[i]} is the weight
  assigned to the \eqn{i}th dissimilarity object and
  \eqn{D_{ijk}}{D[ijk]} is the dissimilarity between \eqn{j} and
  \eqn{k} for the \eqn{i}th dissimilarity object.
}
\value{
  \code{fuse} returns an object of class \code{"dist"} with the
  attribute \code{"method"} set to \code{"fuse"}.

  This is the case even if the supplied objects are full dissimilarity
  matrices. If you want a full dissimilarity object, use
  \code{\link{as.matrix.dist}} on the returned object.

  The returned object contains an extra attribute \code{"weights"},
  which records the weights used in the fusing.
}
\references{
  Melssen W., Wehrens R. and Buydens L. (2006) Supervised Kohonen
  networks for classification problems. \emph{Chemometrics and
  intelligent laboratory systems} \strong{83}, 99--113.
}
\author{Gavin L. Simpson}
\seealso{
  \code{\link{dist}}, \code{\link[vegan]{vegdist}},
  \code{\link{distance}}.
}
\examples{
train1 <- data.frame(matrix(abs(runif(100)), ncol = 10))
train2 <- data.frame(matrix(sample(c(0,1), 100, replace = TRUE),
                     ncol = 10))
rownames(train1) <- rownames(train2) <- LETTERS[1:10]
colnames(train1) <- colnames(train2) <- as.character(1:10)

d1 <- vegdist(train1, method = "bray")
d2 <- vegdist(train2, method = "jaccard")

dd <- fuse(d1, d2, weights = c(0.6, 0.4))
dd
str(dd)
}
\keyword{methods}
\keyword{multivariate}
