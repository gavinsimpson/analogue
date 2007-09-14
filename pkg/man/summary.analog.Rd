\name{summary.analog}
\alias{summary.analog}
\alias{print.summary.analog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarising analogue matching results}
\description{
  \code{\link{summary}} method for class \code{"analog"}.
}
\usage{
\method{summary}{analog}(object, display = c("dist", "names", "quantiles"),
        k = 10, probs = c(0.01, 0.02, 0.05, 0.1, 0.2), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"analog"}, usually as a result
    of a call to \code{\link{analog}}.}
  \item{display}{character; one or more of the listed
    choices. Determines which aspects of the \code{\link{analog}}
    results are summarised.}
  \item{k}{number of analogues to use. If missing, \code{k} is chosen
    automatically as the \code{k} that achieves lowest RMSE.}
  \item{probs}{numeric; giving the probabilities of the distribution to
    return quantiles for. See \code{\link{quantile}}.}
  \item{\dots}{arguments passed to or from other methods.}
}
\value{
  A list with one or more of the components listed below. Attributes
  \code{"method"}, \code{"train"}, \code{"call"} and \code{"k"} contain
  the dissimilarity coefficient used, whether the training set
  dissimilarities were saved, the matched function call and the number
  of close analogues to return respectively.
  
  \item{dists}{a matrix of dissimilarities between training set samples
    and fossil samples. The number of rows is given by argument
    \code{k}. There is a column for each fossil sample.}
  \item{names}{a matrix of names of samples from the training set that
    are analogues for each fossil sample. The number of rows is given by
    argument \code{k}. There is a column for each fossil sample.}
  \item{quantiles}{numeric; the quantiles of the distribution of the
    pairwise dissimilarities for the training set for probabilities
    \code{prob}.}
}
\author{Gavin L. Simpson}
\seealso{\code{\link{analog}}.}
\examples{
\dontrun{
## continue the RLGH example from ?join
example(join)

## analog matching between SWAP and RLGH core
swap.analog <- analog(swapdiat, rlgh, method = "chord")
swap.analog
summary(swap.analog)
}
}
\keyword{methods}% at least one, from doc/KEYWORDS
