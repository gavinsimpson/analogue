\name{minDC}
\alias{minDC}
\alias{minDC.default}
\alias{minDC.predict.mat}
\alias{minDC.analog}
\alias{minDC.wa}
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
minDC(x, \dots)

\method{minDC}{default}(x, \dots)

\method{minDC}{predict.mat}(x, \dots)

\method{minDC}{analog}(x, probs = c(0.01, 0.02, 0.05, 0.1), \dots)

\method{minDC}{wa}(x, y,
      method = c("euclidean", "SQeuclidean", "chord", "SQchord",
                 "bray", "chi.square", "SQchi.square", "information",
                 "chi.distance", "manhattan", "kendall", "gower",
                 "alt.gower", "mixed"),
      percent = FALSE, probs = c(0.01, 0.025, 0.05, 0.1), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"predict.mat"}, \code{"analog"} or
    a object with a component named \code{"minDC"}.}
  \item{probs}{numeric; vector of probabilities with values in [0,1].}
  \item{y}{an optional matrix-like object containing fossil samples for
    which the minimum dissimilarities to training samples are to be
    calculated.}
  \item{method}{character; which choice of dissimilarity coefficient to
    use. One of the listed options. See \code{\link{distance}}.}
  \item{percent}{logical; Are the data percentages? If \code{TRUE},
    the data (\code{x} and \code{y}) will be divided by 100 to convert
    them to the proportions expected by \code{\link{distance}}.}
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

## fit the MAT model using the squared chord distance measure
ik.mat <- mat(ImbrieKipp, SumSST, method = "SQchord")
ik.mat

## reconstruct for the V12-122 core data
v12.mat <- predict(ik.mat, V12.122)

## extract the minimum DC values
v12.mdc <- minDC(v12.mat)
v12.mdc

## draw a plot of minimum DC by time
plot(v12.mdc, use.labels = TRUE, xlab = "Depth (cm.)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{manip}
\keyword{methods}
