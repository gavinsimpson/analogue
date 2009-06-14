\name{stdError}
\alias{stdError}
\alias{stdError.mat}
\alias{stdError.predict.mat}
\title{Standard error of MAT fitted and predicted values}
\description{
  Computes the weighted standard deviation of the environment for the
  \emph{k}-closest analogues for each sample. This measure is proposed
  as a measure of reconstruction uncertainty for MAT models.
}
\usage{
stdError(object, ...)

\method{stdError}{mat}(object, k, ...)

\method{stdError}{predict.mat}(object, k, ...)
}
\arguments{
  \item{object}{Object for which the uncertainty measure is to be
    computed. Currently methods for \code{\link{mat}} and
    \code{\link{predict.mat}}.}
  \item{k}{numeric; how many analogues to take? If missing, the default,
    \code{k} is chosen using \code{\link{getK}}.}
  \item{\dots}{Additional arguments passed to other methods. Currently
    not used.}
}
%\details{
%  TODO
%}
\value{
  A named numeric vector of weighted standard deviations of the
  environment for the \emph{k} closest analogues used to compute the MAT
  predicted values.

  The returned vector has attributes \code{"k"} and \code{"auto"},
  indicating the number of analogues used and whether this was
  determined from \code{object} or supplied by the user.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Gavin L. Simpson}
%\note{ ~~further notes~~ 
%
%}
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

## standard errors
stdError(ik.mat)

## reconstruct for the V12-122 core data
coreV12.mat <- predict(ik.mat, V12.122, k = 3)
## standard errors
stdError(coreV12.mat)
}
\keyword{methods}
\keyword{univar}