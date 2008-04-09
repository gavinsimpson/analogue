\name{performance}
\alias{performance}
\alias{print.performance}
\alias{performance.wa}
\alias{performance.predict.wa}
\alias{performance.bootstrap.wa}
\title{Transfer function model performance statistics}
\description{
  A simple extractor function to access the model performance statistics
  of transfer function models.
}
\usage{
performance(object, ...)

\method{performance}{wa}(object, ...)

\method{performance}{predict.wa}(object, ...)

\method{performance}{bootstrap.wa}(object, ...)
}
\arguments{
  \item{object}{A transfer function object.}
  \item{\dots}{Arguments passed to other methods. Currently ignored.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
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
## continue the example from ?wa
example(wa)

## the model performance statistics
performance(mod)
}
\keyword{methods}
\keyword{utilities}
