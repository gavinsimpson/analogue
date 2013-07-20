\name{caterpillarPlot}
\alias{caterpillarPlot}
\alias{caterpillarPlot.default}
\alias{caterpillarPlot.data.frame}
\alias{caterpillarPlot.wa}
\alias{caterpillar}

\title{
  Caterpillar plot of species' WA optima and tolerance range.
}

\description{
  Draws a caterpillar plot of the weighted average optima and tolerance
  range for each of the species in a data set.
}

\usage{
\method{caterpillarPlot}{default}(x, tol, mult = 1, decreasing = TRUE,
                labels, xlab = NULL, pch = 21, bg = "white",
                col = "black", lcol = col, lwd = 2, frame.plot = FALSE, ...)

\method{caterpillarPlot}{data.frame}(x, env, useN2 = TRUE, xlab, ...)

\method{caterpillarPlot}{wa}(x, type = c("observed","model"), ...)
}


\arguments{
  \item{x}{For the \code{default} method, a numeric vector of species
    optima. For the \code{data.frame} method a species data matrix or data
    frame. For the \code{wa} method an object of class \code{"wa"}.}
  \item{tol}{numeric; vector of species tolerances.}
  \item{env}{numeric; variable for which optima and tolerances are
    required.}
  \item{useN2}{logical; should Hill's N2 values be used to produce
    un-biased tolerances?}
  \item{decreasing}{logical; should the sort order of the species be
    increasing or decreasing?}
  \item{mult}{numeric; multiplication factor for species' tolerances.}
  \item{labels}{character; vector of labels for the species names with
    which to annotate the y-axis. If missing, \code{names(x)} is used.}
  \item{xlab}{character; the x-axis label. If \code{NULL}, the default,
    a description of \code{env} is used.}
  \item{pch, bg, col}{The plotting character to use and its background and
    foreground colour. See \code{\link{par}}.}
  \item{lcol, lwd}{The colour and line width to use for the tolerance
    range.}
  \item{type}{character; \code{"observed"} uses the actual tolerances
    observed from the data. \code{"model"} uses the tolerances used in
    the WA model where very small tolerances have been reset for some
    definition of small.}
  \item{frame.plot}{logical; should a box be drawn round the plot?}
  \item{\dots}{Additional graphical arguments to be passed on to
    plotting functions.}
}

\details{
  The function may also be called using the short form name
  \code{caterpillar}:
  \preformatted{
caterpillar(x, ...)
  }
}

\value{
  The function results in a plot on the currently active device. A data
  frame with components \code{Optima} and \code{Tolerance} is returned
  invisibly.
}

%\references{}

\author{
  Gavin L. Simpson
}

%\note{}

\seealso{
  For the underlying computations \code{\link{optima}} and
  \code{\link{tolerance}}.
}

\examples{
data(ImbrieKipp)
data(SumSST)

## default plot
caterpillar(ImbrieKipp, SumSST)

## customisation
opttol <-
    caterpillar(ImbrieKipp, SumSST, col = "red2",
                bg = "yellow", lcol = "blue",
                xlab = expression(Summer ~ Sea ~ Surface ~
                                 Temperature~(degree*C)))

## invisibly returns the optima and tolerances
head(opttol)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
