\name{residLen}
\alias{residLen}
\alias{print.residLen}
\alias{fittedY}
\alias{sqrlLinear}
\alias{sqrlUnimodal}
\title{Squared residual length diagnostics}
\description{
  The squared residual length between the fitted values of a constrained
  ordination and the original species data is one diagnostic for
  transfer function models.
}
\usage{
residLen(X, env, passive, method = c("cca","rda"))

fittedY(ord, newdata, colsum)

sqrlLinear(Y, fitted)

sqrlUnimodal(Y, colsum, fitted)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{data frame; the training set species data.}
  \item{env}{vector; the training set environmental data.}
  \item{passive}{data frame; the passive samples species data.}
  \item{method}{the ordination technique to use. One of \code{"rda"} or
    \code{"cca"}, with the latter the default.}
  \item{ord}{an ordination object, the result of a call to
    \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.}
  \item{newdata}{Species data matrix for passive samples. Must have same
    columns as data used to fit \code{ord}.}
  \item{colsum}{column (species) sums for training set data used to fit
    \code{ord}.}
  \item{Y}{Original species data matrix, the response for which squared
    residual lengths are to be computed.}
  \item{fitted}{The fitted values of the response derived from the
    constrained ordination model.}
}
\details{
  The squared residual lengths are computed for the training set samples
  and the passive samples separately. Passive samples that are poorly
  fitted in the transfer function model will have large squared residual
  distances between the observed species data and the fitted values from
  the constrained ordination.

  \code{residLen} is the main user-interface function and can be called
  with either the training data and passive samples.

  \code{fittedY} returns the fitted approximation of the passive sample
  response data (i.e. species data). \code{sqrlLinear} and
  \code{sqrlUnimodal} return the squared residual distances between the
  observed species data and the fitted values from the constrained
  ordination model.
}
\value{
  \code{fittedY} returns a matrix of fitted species abundances for
  passive samples.

  \code{sqrlLinear} and \code{sqrlUnimodal} return a vector of
  residual distances.
  
  \code{residLen} returns an object of class \code{"residLen"} with the
  attribute \code{"method"} set to \code{"method"}. This object is a
  list with the following components:
  
  \item{train, passive}{The squared residual lengths for the training
    set and the passive samples.}
  \item{ordination}{The fitted ordination.}
  \item{call}{The matched call.}
}
\references{
  Ter Braak C.J.F. and Smilauer P. (2002) CANOCO Reference manual and
  CanoDraw for Windows User's guide: Software for Canonical Ordination
  (version 4.5). Microcomputer Power (Ithaca, NY, USA), 500pp.
}
\author{Gavin L. Simpson}
\seealso{
  \code{\link[vegan]{cca}} and \code{\link[vegan]{predict.cca}} for some
  of the underlying computations.
}
\examples{
data(swapdiat, swappH, rlgh)

## squared residual lengths for RLGH
rlens <- residLen(swapdiat, swappH, rlgh)
rlens

## as before but using linear RDA
residLen(swapdiat, swappH, rlgh, method = "rda")
}
\keyword{methods}
\keyword{multivariate}
