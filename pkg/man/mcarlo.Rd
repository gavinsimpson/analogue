\name{mcarlo}
\alias{mcarlo}
\alias{mcarlo.default}
\alias{mcarlo.mat}
\alias{mcarlo.analog}
\alias{print.mcarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Monte Carlo simulation of dissimilarities}
\description{
  Permutations and Monte Carlo simulations to define critical values for
  dissimilarity coefficients for use in MAT reconstructions.
}
\usage{
mcarlo(object, ...)

\method{mcarlo}{default}(object, nsamp = 10000,
       type = c("paired", "complete", "bootstrap", "permuted"),
       replace = FALSE, 
       method = c("euclidean", "SQeuclidean", "chord", "SQchord",
                  "bray", "chi.square", "SQchi.square",
                  "information", "chi.distance", "manhattan",
                  "kendall", "gower", "alt.gower", "mixed"),
       is.dcmat = FALSE, diag = FALSE, ...)

\method{mcarlo}{mat}(object, nsamp = 10000,
       type = c("paired", "complete", "bootstrap", "permuted"),
       replace = FALSE, diag = FALSE, ...)

\method{mcarlo}{analog}(object, nsamp = 10000,
       type = c("paired", "complete", "bootstrap", "permuted"),
       replace = FALSE, diag = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object. Currently only object's of class
    \code{"mat"}, \code{"analog"} or matrix-like object of species data
    allowed.}
  \item{nsamp}{numeric; number of permutations or simulations to draw.}
  \item{type}{character; the type of permutation or simulation to
    perform. See Details, below.}
  \item{replace}{logical; should sampling be done with replacement?}
  \item{method}{character; for raw species matrices, the dissimilarity
    coefficient to use. This is predefined when fitting a MAT model with
    \code{\link{mat}} or analogue matching via \code{\link{analogue}}
    and is ignored in the \code{"mcarlo"} methods for classes
    \code{"mat"} and \code{"analog"}.}
  \item{is.dcmat}{logical; is \code{"object"} a dissimilarity
    matrix. Not meant for general use; used internally by \code{"mat"}
    and \code{"analogue"} methods to instruct the \code{"default"}
    method that \code{"object"} is already a dissimilarity matrix, so
    there is no need to recalculate.}
  \item{diag}{logical; should the dissimilarities include the diagonal
    (zero) values of the dissimilarity matrix. See Details.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  Only \code{"type"} \code{"paired"} and \code{"bootstrap"} are
  currently implemented.

  \code{\link{distance}} produces square, symmetric
  dissimilarity matrices for training sets. The upper triangle of these
  matrices is a duplicate of the lower triangle, and as such is
  redundant. \code{mcarlo} works on the lower triangle of these
  dissimilarity matrices, representing all pairwise dissimilarity values
  for training set samples. The default is \strong{not} to include the
  diagonal (zero) values of the dissimilarity matrix. If you feel that
  these diagonal (zero) values are part of the population of
  dissimilarities then use \code{"diag = TRUE"} to include them in the
  permutations.
}
\value{
  A vector of simulated dissimilarities of length \code{"nsamp"}. The
  \code{"method"} used is stored in attribute \code{"method"}.
}
\references{
  
  Sawada, M., Viau, A.E., Vettoretti, G., Peltier, W.R. and Gajewski,
  K. (2004) Comparison of North-American pollen-based temperature and
  global lake-status with CCCma AGCM2 output at 6 ka. \emph{Quaternary
    Science Reviews} \strong{23}, 87--108.

}
\author{Gavin L. Simpson}
\note{
  The performance of these permutation and simulation techniques still
  needs to be studied. This function is provided for pedagogic
  reasons. Although recommended by Sawada et al (2004), sampling with
  replacement (\code{"replace = TRUE"}) and including diagonal (zero)
  values (\code{"diag = TRUE"}) simulates too many zero distances. This
  is because the same training set sample can, on occasion be drawn
  twice leading to a zero distance. It is impossible to find in nature
  two samples that will be perfectly similar, and as such sampling
  \strong{with} replacement \strong{and} \code{"diag = TRUE"} seems
  undesirable at best.
}
\seealso{\code{\link{mat}} for fitting MAT models and
  \code{\link{analog}} for analogue matching.
  \code{\link{roc}} as an alternative method for determining critical
  values for dissimilarity measures when one has grouped data.
  
  \code{\link{plot.mcarlo}} provides a plotting method to visualise the
  distribution of simulated dissimilarities.
}
\examples{
## continue from example(join)
example(join)

## perform the modified method of Sawada (2004) - paired sampling,
## with replacement
swap.mcarlo <- mcarlo(swapdiat, method = "SQchord", nsamp = 1000,
                      type = "paired", replace = FALSE)
swap.mcarlo

## plot the simulated distribution
par(mfrow = c(2,1))
plot(swap.mcarlo)
par(mfrow = c(1,1))
}
\keyword{multivariate}
%\keyword{ ~kwd2 } % more keywords TO DO
