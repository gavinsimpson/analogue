\name{dissimilarities}
\alias{dissimilarities}
\alias{dissimilarities.analog}
\alias{dissimilarities.mat}
\alias{dissim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract dissimilarity coefficients from models }
\description{
  Extracts a vector of dissimilarity coefficients from an object for
  further analysis.
}
\usage{

dissimilarities(object, \dots)
dissim(object, \dots)

\method{dissimilarities}{analog}(object, which = c("train", "analogs"),
                 \dots)

\method{dissimilarities}{mat}(object, \dots)
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
## Imbrie and Kipp example
## load the example data
data(ImbrieKipp)
data(SumSST)
data(V12.122)

## merge training and test set on columns
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- dat[[1]] / 100
V12.122 <- dat[[2]] / 100

## analog matching between SWAPImbrie & Kipp and V12.122 core
ik.analog <- analog(ImbrieKipp, V12.122, method = "chord")
ik.analog
summary(ik.analog)

## compare training set dissimilarities with normals
## and derive cut-offs
ik.dij <- dissim(ik.analog)
plot(ik.dij)

}
\keyword{methods}
\keyword{multivariate}
