\name{performance}
\alias{performance}
\alias{print.performance}
\alias{performance.wa}
\alias{performance.predict.wa}
\alias{performance.bootstrap.wa}
\alias{performance.crossval}
\title{Transfer function model performance statistics}
\description{
  A simple extractor function to access the model performance statistics
  of transfer function models.
}
\usage{
performance(object, ...)
}
\arguments{
  \item{object}{A transfer function object.}
  \item{\dots}{Arguments passed to other methods. Currently ignored.}
}
\details{
  \code{performance} is a generic function for use with a number of
  fitted models objects in \pkg{analogue}. The available methods are:

  \describe{
    \item{\code{\link{wa}} }{Weighted Averaging Models.}
    \item{\code{\link{predict.wa}} }{Predictions from a Weighted Average
      Model.}
    \item{\code{\link{pcr}} }{Principal Component Regression models.}
    \item{\code{\link{bootstrap.wa}} }{Bootstrapped Weighted Averaging
      Models.}
    \item{\code{\link{crossval}} }{Cross-validated models fitted via
      \code{\link{crossval}}.}
  }
}
\value{
  A named vector containing the extracted model performance statistics.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Gavin L. Simpson }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{wa}}, \code{\link{predict.wa}},
  \code{\link{bootstrap.wa}}.}
\examples{
data(ImbrieKipp)
data(SumSST)

## fit the WA model
mod <- wa(SumSST ~., data = ImbrieKipp)
mod

## the model performance statistics
performance(mod)
}
\keyword{methods}
\keyword{utilities}
