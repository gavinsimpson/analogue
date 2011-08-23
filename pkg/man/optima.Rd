\name{optima}
\alias{optima}
\alias{optima.default}
\alias{print.optima}
\alias{print.tolerance}
\alias{tolerance.default}
\alias{as.data.frame.optima}
\alias{as.data.frame.tolerance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted averaging optima and tolerance ranges}
\description{
  Computes weighted average optima and tolerance ranges from species
  abundances and values of the environment.
}
\usage{
optima(x, ...)

\method{optima}{default}(x, env, ...)

\method{tolerance}{default}(x, env, useN2 = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Species data matrix or data frame.}
  \item{env}{Numeric; variable for which optima or tolerances are
    required.}
  \item{useN2}{logical; should Hill's N2 values be used to produce
    un-biased tolerances?}
  \item{\dots}{Arguments passed to other methods.}
}
%\details{
%
%}
\value{
  Both functions return a named vector containing the WA optima or
  tolerances for the environmental gradient specified by \code{env}.
}
\references{TO DO}
\author{Gavin L. Simpson}
\note{
  Objects of class \code{"optima"} or \code{"tolerance"} can be coerced
  to data frames using methods for \code{link{as.data.frame}}.
}
\seealso{\code{\link{wa}}}
\examples{
## Load the Imbrie & Kipp data and
## summer sea-surface temperatures
data(ImbrieKipp)
data(SumSST)

## WA optima
(opt <- optima(ImbrieKipp, SumSST))

## WA tolerances
(tol <- tolerance(ImbrieKipp, SumSST, useN2 = TRUE))

## convert to data frame
as.data.frame(opt)
as.data.frame(tol)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
