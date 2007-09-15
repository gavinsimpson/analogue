\name{summary.cma}
\alias{summary.cma}
\alias{print.summary.cma}
\title{Summarise the extraction of close modern analogues}
\description{
  \code{\link{summary}} method for class \code{"cma"}.
}
\usage{
\method{summary}{cma}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{"cma"}, usually the
    result of a call to \code{\link{cma}}.}
  \item{\dots}{arguments passed to or from other methods.}
}
\value{
  An object of class \code{"summary.cma"} with the components of an
  object of class \code{\link{cma}}, plus:  
  \item{distances}{a matrix of distances/dissimilarities. Individual
    columns contain the ordered close modern analogues for individual
    fossil samples. Rows of this matrix refer to the
    \eqn{k^{\mathrm{th}}}{k-th} closest analogue for each fossil
      sample. See notes below.}
  \item{samples }{a matrix of sample names from the reference set that
    are close modern analogues for a fossil sample. Individual
    columns contain the ordered close modern analogues for individual 
    fossil samples. Rows of this matrix refer to the
    \eqn{k^{\mathrm{th}}}{k-th} closest analogue for each fossil sample.
      See notes below.}
}
\author{Gavin L. Simpson}
\note{
  Currently, only objects of class \code{\link{analog}} are supported.
  The number of rows in the returned matrices is equal to the maximum
  number of close modern analogues identified for an individual fossil
  sample. If no close modern analogues exist for an individual fossil
  sample are identified, then the relevant column in \code{"distances"}
  will contain all missing values and in \code{"samples"} the string
  \code{"none"}. Rows of individual columns will be padded with missing
  values if the number of close modern analogues for that sample is less
  than the maximum number of close modern analogues identified for a
  single sample.
}
\seealso{\code{\link{cma}}}
\examples{
\dontrun{
## continue the RLGH example from ?join
example(join)

## analog matching between SWAP and RLGH core
swap.analog <- analog(swapdiat, rlgh, method = "chord")
swap.analog
summary(swap.analog)

## close modern analogues
swap.cma <- cma(swap.analog, cutoff = 0.6)
swap.cma
summary(swap.cma)
}
}
\keyword{methods}
