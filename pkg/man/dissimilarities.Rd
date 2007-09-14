\name{dissimilarities}
\alias{dissimilarities}
\alias{dissimilarities.analog}
\alias{dissim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract dissimilarity coefficients from models }
\description{
  Extracts a vector of dissimilarity coefficients from an object for
  further analysis.
}
\usage{

dissimilarities(object, ...)
dissim(object, ...)

\method{dissimilarities}{analog}(object, which = c("train", "analogs"),
                 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object from which the dissimilarity values are to
    be extracted. Currently only for objects of class \code{"analog"}.}
  \item{which}{character; which set of dissimilarities should be
    extracted. One of \code{"train"} or \code{"analogs"}.}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  The function can be called using the much shorter name
  \code{"dissim"}.
}
\value{
  A vector of dissimilarities.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{analog}}, \code{\link{plot.dissimilarities}} }
\examples{
## continue the RLGH example from ?join
example(join)

## analog matching between SWAP and RLGH core
swap.analog <- analog(swapdiat, rlgh, method = "chord")
swap.analog
summary(swap.analog)

## compare training set dissimilarities with normals
## and derive cut-offs
swap.dissim <- dissim(swap.analog)
plot(swap.dissim)

}
\keyword{methods}
\keyword{multivariate}
