\name{residLen}
\alias{residLen}
\alias{residLen.default}
\alias{print.residLen}
\title{Squared residual length diagnostics}
\description{
  The squared residual length between the fitted values of a constrained
  ordination and the original species data is one diagnostic for
  transfer function models.
}
\usage{
residLen(train, ...)

\method{residLen}{default}(train, env, passive,
         method = c("cca","rda"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{data frame; the training set species data.}
  \item{env}{vector; the training set environmental data.}
  \item{passive}{data frame; the passive samples species data.}
  \item{method}{the ordination technique to use. Only \code{"cca"} is
    currently supported.}
  \item{\dots}{additional arguments passed to methods.}
}
\details{
  The squared residual lengths are computed for the training set samples
  and the passive samples separately. Passive samples that are poorly
  fitted in the transfer function model will have large squared residual
  distances between the observed species data and the fitted values from
  the constrained ordination.
}
\value{
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
}
\keyword{methods}
\keyword{multivariate}
