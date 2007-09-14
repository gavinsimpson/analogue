\name{screeplot}
\alias{screeplot.mat}
\alias{screeplot.bootstrap.mat}
%\alias{Screeplot.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Screeplots of model results}
\description{
  Draws screeplots of performance statistics for models of varying complexity.
}
\usage{

\method{screeplot}{mat}(x, k, restrict = 20,
          display = c("rmsep", "avg.bias",
                      "max.bias", "r.squared"),
          weighted = FALSE,  col = "red", xlab = NULL,
          ylab = NULL, main = NULL, sub = NULL, \dots)

\method{screeplot}{bootstrap.mat}(x, k, restrict = 20,
          display = c("rmsep","avg.bias","max.bias",
                      "r.squared"),
          legend = TRUE, loc.legend = "topright",
          col = c("red", "blue"),
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          \dots,
          lty = c("solid","dashed"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{\link{mat}} and \code{bootstrap.mat}.}
  \item{k}{number of analogues to use. If missing 'k' is chosen
    automatically as the 'k' that achieves lowest RMSE.} 
  \item{restrict}{logical; restrict comparison of k-closest model to k
    \eqn{<=}{<=} \code{restrict}.}
  \item{display}{which aspect of \code{x} to plot? Partial match.}
  \item{weighted}{logical; should the analysis use weighted mean of env
    data of analogues as fitted/estimated values?}
  \item{xlab, ylab}{x- and y-axis labels respectively.}
  \item{main, sub}{main and subtitle for the plot.}
  \item{legend}{logical; should a legend be displayed on the figure?}
  \item{loc.legend}{character; a keyword for the location of the
    legend. See \code{\link{legend}} for details of allowed keywords.}
  \item{col}{Colours for lines drawn on the screeplot. Method for class
    \code{"bootstrap.mat"} takes a vector of two colours.}
  \item{lty}{vector detailing the line type to use in drawing the
    screeplot of the apparent and bootstrap statistics,
    respectively. Code currently assumes that \code{length(lty)} is 2.}
  \item{\dots}{arguments passed to other graphics functions.}
}
\details{
  Screeplots are often used to graphically show the results of
  cross-validation or other estimate of model performance across a range
  of model complexity.

  Four measures of model performance are currently available: i) root
  mean square error of prediction (RMSEP); ii) average bias --- the
  mean of the model residuals; iii) maximum bias --- the maximum average
  bias calculated for each of \emph{n} sections of the gradient of the
  environmental variable; and v) model \eqn{R^2}.

  For the maximum bias statistic, the response (environmental) gradient
  is split into \emph{n} = 10 sections.

  For the \code{\link{bootstrap}} method, apparent and bootstrap
  versions of these statistics are available and plotted.
}
\author{Gavin Simpson}
\note{
  Currently only models of class \code{\link{mat}} and
  \code{\link{bootstrap.mat}} are supported.
}
\seealso{\code{\link{screeplot}} }
\examples{
## continue the example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

##
screeplot(swap.mat)
}
\keyword{hplot}
\keyword{methods}
