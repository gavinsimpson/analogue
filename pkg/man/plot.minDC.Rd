\name{plot.minDC}
\alias{plot.minDC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of minimum dissimilarity per sample}
\description{
  Minimum dissimilarity is a useful indicator of reliability of
  reconstructions performed via MAT and other methods, and for analogue
  matching. Minimum dissimilarity for a sample is the smallest
  dissimilarity between it and the training set samples.
}
\usage{
\method{plot}{minDC}(x, depths, use.labels = FALSE,
           quantiles = TRUE, rev.x = TRUE, type = "l",
           xlim, ylim, xlab = "", ylab = "Dissimilarity",
           main = "", sub = NULL,
           col.quantile = "red", lty.quantile = "dotted",
           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"minDC"}.}
  \item{depths}{numeric; a vector of depths for which predicted values
    exist or will be generated. Can be missing, in which case,
    \strong{if} \code{use.labels = TRUE}, the function will attempt to
    derive suitable values for you. See Details below.}
  \item{use.labels}{logical; should \code{reconPlot} attempt to derive
    values for argument \code{depths} from the names of the predicted
    values? Only use if \code{depths} is missing. See Details below.}
  \item{quantiles}{logical; should the probability quantiles be drawn on
  the plot?}
  \item{rev.x}{logical; should the depth/age axis be reversed (drawn
    from high to low)?}
  \item{type}{type of line drawn. See \code{\link{par}} and argument
    \code{"type"}.}
  \item{xlab, ylab}{character; the x- and y-axis labels respectively.}
  \item{main, sub}{character; main title and subtitle for the plot.}
  \item{xlim, ylim}{numeric, length 2; the x- and y-limits for the
    plotted axes. If not provided, the function will calculate
    appropriate values to cover the range of plotted values and any
    quantile lines (if requested via \code{"quantiles"}.}
  \item{col.quantile}{colour in which to draw the quantile lines.}
  \item{lty.quantile}{line type in which to draw the quantile lines.}
  \item{\dots}{arguments to be passed to methods, such as graphical
          parameters (see \code{\link{par}}). Currently ignored.}
}
\details{
  Conventionally, these plots are drawn on a depth or an age
  scale. Argument \code{depths} is used to provide the depth or age
  axis, against which the predicted values are plotted.

  If \code{depths} is not provided, then the function will try to
  derive the appropriate values from the labels of the predictions if
  \code{use.labels = TRUE}. You must provide \code{depths} or set
  \code{use.labels = TRUE} otherwise an error will result. The derived
  labels will be coerced to numerics. If your labels are coercible, then
  you'll either get nonsense on the plot or an error from R. If so,
  provide suitable values for \code{depths}.
}
\value{
  A plot on the currently active device.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Gavin L. Simpson}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{minDC}}}
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
\keyword{hplot}% __ONLY ONE__ keyword per line
\keyword{methods}
