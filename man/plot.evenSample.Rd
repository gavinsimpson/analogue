\name{plot.evenSample}
\alias{plot.evenSample}

\title{Plot distribution of samples along gradient}
\description{
  \code{\link{plot}} method for objects produced by
  \code{\link{evenSample}}. Draws a histogram of the number of samples
  per gradient segment.
}
\usage{
\method{plot}{evenSample}(x, add = FALSE, xlim = NULL, ylim = NULL, col = "grey",
     border = "grey", lty = NULL, ylab, xlab, main = NULL, sub = NULL,
     ann = TRUE, axes = TRUE, \dots)
}

\arguments{
  \item{x}{an object of class \code{"evenSample"}, resulting from a call
    to \code{\link{evenSample}}}
  \item{add}{logical; should the histogram of counts be added to an
    existing plot?}
  \item{xlim, ylim}{numeric; user-specified axis limits. If not suplied
    suitable defaults are used}
  \item{col, border}{colours for the fill and border of the histogram
    bars respectively.}
  \item{lty}{the line type with which to draw the histogram bars}
  \item{ylab, xlab, main, sub}{character strings used to label the plot}
  \item{ann}{logical; should the default annotations be added to the
    plot. This relates to the plot main and sub titles and the x and y
    axis labels.}
  \item{axes}{logical; should plot axes be drawn?}
  \item{\dots}{additional arguments passed to/from other methods.}
}

\value{
  A plot is draw is drawn on the currently active device.
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Gavin L. Simpson
}

\seealso{
  \code{\link{evenSample}} for the function used to create objects used
  by this plot method.
}
\keyword{ hplot }
\keyword{ method }
