\name{predict.wa}
\alias{predict.wa}
\alias{print.predict.wa}
\title{Predict from a weighted average model}
\description{
  Model predictions and cross-validation predictions for weighted
  averaging transfer function models.
}
\usage{
\method{predict}{wa}(object, newdata,
        CV = c("none", "LOO", "bootstrap"),
        verbose = FALSE, n.boot = 100, ...)
}
\arguments{
  \item{object}{an object of class \code{"wa"}, usually the result of a
    call to \code{\link{wa}}}
  \item{newdata}{An optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{CV}{Should cross-validation be performed? Leave-one-out
    (\code{"LOO"}) and bootstrap (\code{"bootstrap"}) CV are currently
    available.}
  \item{verbose}{Should CV progress be printed to the console?}
  \item{n.boot}{The number of bootstrap smaples to draw.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Not all CV methods produce the same output. \code{CV = "bootstrap"} is
  currently the only method that produces sample specific errors.
}
\value{
  An object of class \code{"predict.wa"}, a list with the following
  components:

  \item{pred}{A list with components \code{pred} and \code{rmsep}
    containing the predicted values and the sample specific errors if
    available.}
  \item{performance}{A list with model performance statistics.}
  \item{model.pred}{A list with components \code{pred} and \code{rmsep}
    containing the predicted values for the training set samples and the
    sample specific errors if available.}
  \item{call}{the matched function call.}
  \item{CV.method}{The CV method used.}
}
\references{
  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C. and ter Braak,
  C.J.F. (1990). Diatoms and pH reconstruction. \emph{Philosophical
    Transactions of the Royal Society of London; Series B}, \strong{327};
  263--278.}
\author{Gavin L. Simpson}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{wa}}.}
\examples{
## continue the example from ?wa
example(wa)

## load the RLGH data
data(rlgh)
rlgh <- rlgh / 100

## Predict pH for the RLGH samples
rlgh.pred <- predict(mod, rlgh, CV = "bootstrap", n.boot = 100)

## draw the fitted reconstruction
reconPlot(rlgh.pred, use.labels = TRUE, display = "bars")

## extract the model performance stats
performance(rlgh.pred)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
