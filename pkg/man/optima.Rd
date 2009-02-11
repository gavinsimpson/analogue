\name{optima}
\alias{optima}
\alias{optima.default}
\alias{print.optima}
\alias{tolerance}
\alias{print.tolerance}
\alias{tolerance.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted averaging optima and tolerance ranges}
\description{
  Computes weighted average optima and tolerance ranges from species
  abundances and values of the environment.
}
\usage{
optima(x, ...)

\method{optima}{default}(x, env, ...)

tolerance(x, ...)

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
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Both functions return a named vector containing the WA optima or
  tolerances for the environmental gradient specified by \code{env}.
}
\references{TO DO}
\author{Gavin L. Simpson}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{wa}}}
\examples{
## Load the Imbrie & Kipp data and
## summer sea-surface temperatures
data(ImbrieKipp)
data(SumSST)

## WA optima
optima(ImbrieKipp, SumSST)

## WA tolerances
tolerance(ImbrieKipp, SumSST, useN2 = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
