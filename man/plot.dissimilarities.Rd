\name{plot.dissimilarities}
\alias{plot.dissimilarities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the distribution of extracted dissimilarities}
\description{
  Produces a plot of the distribution of the extracted dissimilarities
  and a reference normal distribution with comparable mean and sd.
}
\usage{
\method{plot}{dissimilarities}(x, prob = 0.05,
     legend = TRUE, n.rnorm = 1e+05, col = "black",
     col.ref = "red", lty = "solid", lty.quant = "dotted",
     xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"dissimilarities"}.}
  \item{prob}{numeric; density probability defining the threshold for
    close modern analogues.}
  \item{legend}{logical; draw a legend on the plotted figure?}
  \item{n.rnorm}{numeric; number of random normal deviates for reference
    line.}
  \item{col, col.ref}{colours for the dissimilarity and reference
    density functions drawn.}
  \item{lty, lty.quant}{line types for the dissimilarity and reference
    density functions drawn.}
  \item{xlab, ylab}{character; x- and y-axis labels.}
  \item{main, sub}{character; main and subtitle for the plot.}
  \item{\dots}{graphical arguments passed to other graphics functions.}
}
%\details{
%  FIXME: TO DO
%}
\value{
  A plot on the currently active device.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{dissimilarities}}}
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
\keyword{hplot}% at least one, from doc/KEYWORDS
