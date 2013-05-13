\name{scores.prcurve}
\alias{scores.prcurve}

\title{\code{\link[vegan]{scores}} method for principal curve objects of
  class \code{"prcurve"}.}
\description{
  A \code{\link[vegan]{scores}} method to extract the position on the
  curve to which each observation projects (\code{display = "curve"}) or
  the coordinates of the curve in the dimensions of the input data
  (\code{display = "dimensions"}).
}
\usage{
\method{scores}{prcurve}(x, display = c("curve", "dimensions"), ...)
}

\arguments{
  \item{x}{an object of class \code{"prcurve"}, usually from a call to
    \code{\link{prcurve}}.}
  \item{display}{character; which type of scores to
    extract. \code{display = "curve"} returns the position along the curve
    onto which each observation projects; this can be used like a PCA axis
    score. \code{display = "dimensions"} returns the coordinates of the
    curve in the dimensions of the original data.}
  \item{\dots}{Arguments passed to other methods. Not used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  If \code{display = "curve"} a 1-column matrix is returned with a row
  for each observation in the input data. If \code{display =
  "dimensions"}, a matrix of coordinates for the principal curve is
  returned. The dimensions of this matrix relate to the dimensions of
  the input data; if there were \eqn{n} samples (rows) and \eqn{m}
  variables {columns} then the matrix returned by \code{scores.prcurve}
  will have \eqn{n} rows and \eqn{m} columns.
}
\author{Gavin L. Simpson}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{prcurve}} for fitting principal curves to data.
}
\examples{
## Load the Abernethy Forest data set
data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit the principal curve using varying complexity of smoothers
## for each species
aber.pc <- prcurve(abernethy2, method = "ca", trace = TRUE,
                   vary = TRUE, penalty = 1.4)

## Extract position on the curve
pos <- scores(aber.pc, display = "curve")
head(pos)

## Extract the coordinates of the curve
coord <- scores(aber.pc, display = "dimensions")
dim(coord)
all.equal(dim(coord), dim(abernethy2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
