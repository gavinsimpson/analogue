\name{stdError}
\alias{stdError}
\alias{stdError.mat}
\alias{stdError.predict.mat}
\title{Standard error of MAT fitted and predicted values}
\description{
  Computes the (weighted) standard deviation of the environment for the
  \emph{k}-closest analogues for each sample. This was proposed as one
  measure of reconstruction uncertainty for MAT models (ter Braak,
  1995).
}
\usage{
stdError(object, ...)

\method{stdError}{mat}(object, k, weighted = FALSE, ...)

\method{stdError}{predict.mat}(object, k, weighted = FALSE, ...)
}
\arguments{
  \item{object}{Object for which the uncertainty measure is to be
    computed. Currently methods for \code{\link{mat}} and
    \code{\link{predict.mat}}.}
  \item{k}{numeric; how many analogues to take? If missing, the default,
    \code{k} is chosen using \code{\link{getK}}.}
  \item{weighted}{logical; use a weighted computation?}
  \item{\dots}{Additional arguments passed to other methods. Currently
    not used.}
}
\details{
  Two types of standard error can be produced depending upon whether the
  mean or weighted mean of \eqn{y} for the \eqn{k} closest analogues is
  used for the MAT predictions. If \code{weighted = FALSE} then the
  usual standard deviation of the response for the \eqn{k} closest
  analogues is returned, whereas for \code{weighted = TRUE} a weighted
  standard deviation is used. The weights are the inverse of the
  dissimilarity between the target observation and each of the \eqn{k}
  closest analogues.
}
\value{
  A named numeric vector of weighted standard deviations of the
  environment for the \emph{k} closest analogues used to compute the MAT
  predicted values.

  The returned vector has attributes \code{"k"} and \code{"auto"},
  indicating the number of analogues used and whether this was
  determined from \code{object} or supplied by the user.
}
\references{
  Simpson, G.L. (2012) Analogue methods in palaeolimnology. In Birks,
  H.J.B, Lotter, A.F. Juggins S., and Smol, J.P. (Eds) \emph{Tracking
  Environmental Change Using Lake Sediments, Volume 5: Data Handling and
  Numerical Techniques}. Springer, Dordrecht. 
  
  ter Braak, C.J.F. (1995) Non-linear methods for multivariate
  statistical calibration and their use in palaeoecology: a comparison
  of inverse (\emph{k}-nearest neighbours, partial least squares, and
  weighted averaging partial least squares) and classical
  approaches. \emph{Chemometrics and Intelligent Laboratory Systems}
  \strong{28}:165--180.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{minDC}}, \code{\link{mat}},
  \code{\link{predict.mat}}.}
\examples{
## Imbrie and Kipp Sea Surface Temperature
data(ImbrieKipp)
data(SumSST)
data(V12.122)

## merge training set and core samples
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- dat[[1]] / 100
ImbrieKippCore <- dat[[2]] / 100

## fit the MAT model using the squared chord distance measure
ik.mat <- mat(ImbrieKipp, SumSST, method = "SQchord")

## standard errors - unweighted
stdError(ik.mat)
## standard errors - weighted version for above
stdError(ik.mat, k = getK(ik.mat), weighted = TRUE)

## standard errors - weighted; note this uses more (7) analogues
## than the above as this model had lowest LOO error
stdError(ik.mat, weighted = TRUE)

## reconstruct for the V12-122 core data
coreV12.mat <- predict(ik.mat, V12.122, k = 3)
## standard errors
stdError(coreV12.mat)
}
\keyword{methods}
\keyword{univar}