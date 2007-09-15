\name{minDC}
\alias{minDC}
\alias{minDC.default}
\alias{minDC.predict.mat}
\alias{minDC.analog}
\alias{print.minDC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract minimum dissimilarities}
\description{
  Minimum dissimilarity is a useful indicator of reliability of
  reconstructions performed via MAT and other methods, and for analogue
  matching. Minimum dissimilarity for a sample is the smallest
  dissimilarity between it and the training set samples.
}
\usage{
minDC(x, ...)

\method{minDC}{default}(x, ...)

\method{minDC}{predict.mat}(x, ...)

\method{minDC}{analog}(x, probs = c(0.01, 0.02, 0.05, 0.1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"predict.mat"}, \code{"analog"} or
    a object with a component named \code{"minDC"}.}
  \item{probs}{numeric; vector of probabilities with values in [0,1].}
  \item{\dots}{other arguments to be passed to other methods. Currently
    ignored.}
}
%\details{
%}
\value{
  \code{minDC} returns an object of class \code{"minDC"}.

  An object of class \code{minDC} is a list with some or all of the
  following components:
  
  \item{minDC }{numeric; vector of minimum dissimilarities.}
  \item{method }{character; the dissimilarity coefficient used.}
  \item{quantiles }{numeric; named vector of probability quantiles for
    distribution of dissimilarities of modern training set.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Gavin L. Simpson}
\note{
The \code{"default"} method of \code{minDC} will attempt to extract the
relevant component of the object in question. This may be useful until a
specific \code{minDC} method is written for a given class.
}
\seealso{\code{\link{predict.mat}}, and \code{\link{plot.minDC}} for a
  plotting method.}
\examples{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

## reconstruct for the RLGH core data
rlgh.mat <- predict(swap.mat, rlgh, k = 10)

## extract the minimum DC values
rlgh.mdc <- minDC(rlgh.mat)
rlgh.mdc

## draw a plot of minimum DC by time
plot(rlgh.mdc, use.labels = TRUE, xlab = "Depth (cm.)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{manip}
\keyword{methods}
