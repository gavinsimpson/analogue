\name{reconPlot}
\alias{reconPlot}
\alias{reconPlot.default}
\alias{reconPlot.predict.mat}
\alias{reconPlot.predict.wa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stratigraphic plots of palaeoenvironmental reconstructions}
\description{
  Draws a palaeoenvironmental reconstruction of predicted environmental
  values for sub-fossil assemblages.
}
\usage{

reconPlot(x, \dots)

\method{reconPlot}{default}(x, depths, errors,
          display.error = c("none", "bars", "lines"),
          rev.x = TRUE,
          col.error = "grey", lty.error = "dashed",
          type = "l",
          xlim, ylim,
          xlab = "", ylab = "", main = "",
          \dots)

\method{reconPlot}{predict.mat}(x, depths, use.labels = FALSE,
          predictions = c("model", "bootstrap"),
          display.error = c("none", "bars", "lines"),
          sample.specific = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An R object.}
  \item{depths}{numeric; a vector of depths for which predicted values
    exist or will be generated. Can be missing, in which case,
    \strong{if} \code{use.labels = TRUE}, the function will attempt to
    derive suitable values for you. See Details below.}
  \item{errors}{numeric; a vector of errors for plotting error bars or
    lines.}
  \item{display.error}{character; hown should error bars be drawn on the
    plot? One of \code{"none"}, \code{"bars"}, or \code{"lines"}. If
    \code{"bars"}, error bars are drawn for each sample. If
    \code{"lines"}, lines are drawn enclosing the region prediction +/-
    RMSEP.}
  \item{rev.x}{logical; should the depth/age axis be reversed (drawn
    from high to low)?}
  \item{col.error, lty.error}{the colour and type of line drawn. See
    \code{\link{par}} and arguments \code{"col"} and \code{"lty"}.}
  \item{type}{type of line drawn. See \code{\link{par}} and argument
    \code{"type"}.}
  \item{xlab, ylab}{character; the x- and y-axis labels respectively.}
  \item{main}{character; main title for the plot.}
  \item{xlim, ylim}{numeric, length 2; the x- and y-limits for the
    plotted axes. If not provided, the function will calculate
    appropriate values to cover the range of plotted values and any
    error bars (if requested via \code{"error.bars"}.}
  \item{use.labels}{logical; should \code{reconPlot} attempt to derive
    values for argument \code{depths} from the names of the predicted
    values? Only use if \code{depths} is missing. See Details below.}
  \item{predictions}{character; one of \code{"model"} or
    \code{"bootstrap"}. Which type of predicted values should be
    plotted? The model predictions (\code{"model"}) or the
    bootstrap-derived predictions (\code{"bootstrap"}).}
  \item{sample.specific}{logical; should sample specific errors be used?
  Only for \code{predictions = "bootstrap"}.}
  \item{\dots}{arguments passed to other graphics functions.}
}
\details{
  Conventionally, these plots are drawn on a depth or an age
  scale. Argument \code{depths} is used to provide the depth or age
  axis, against which the predicted values are plotted.

  If \code{depths} is not provided, then the function will try to
  derive the appropriate values from the labels of the predictions if
  \code{use.labels = TRUE}. You must provide \code{depths} or set
  \code{use.labels = TRUE} otherwise an error will result. The derived
  labels will be coerced to numerics. If your labels are not coercible,
  then you'll either get nonsense on the plot or an error from R. If so,
  provide suitable values for \code{depths}.
}
\value{
  A plot on the currently active device.
}
\author{Gavin L. Simpson }
\seealso{\code{\link{mat}}, and \code{\link{predict.mat}} }
\examples{
## Continue example from ?join
example(join)

## Fit a MAT model
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

## Reconstruct pH for the RLGH core
rlgh.pH <- predict(swap.mat, rlgh)

## draw the reconstruction
reconPlot(rlgh.pH, use.labels = TRUE, display.error = "bars",
          xlab = "Depth", ylab = "pH")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
