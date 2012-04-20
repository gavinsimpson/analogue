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
        CV = c("none", "LOO", "bootstrap", "nfold"),
        verbose = FALSE, n.boot = 100, nfold = 5, ...)
}
\arguments{
  \item{object}{an object of class \code{"wa"}, usually the result of a
    call to \code{\link{wa}}}
  \item{newdata}{An optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{CV}{Should cross-validation be performed? Leave-one-out
    (\code{"LOO"}), bootstrap (\code{"bootstrap"}) and \eqn{k}-fold
    (\code{"nfold"}) CV are currently available.}
  \item{verbose}{Should CV progress be printed to the console?}
  \item{n.boot}{The number of bootstrap samples or \eqn{k}-fold steps.}
  \item{nfold}{Number of subsets in \eqn{k}-fold CV.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Not all CV methods produce the same output. \code{CV = "bootstrap"} and
  \code{CV = "nfold"} produce sample specific errors.
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
\author{Gavin L. Simpson and Jari Oksanen (\eqn{k}-fold CV)}
\seealso{\code{\link{wa}}, \code{\link{predict.mat}},
    \code{\link{performance}}, \code{\link{reconPlot}}.}
\examples{
## Imbrie and Kipp
data(ImbrieKipp)
data(SumSST)
ik.wa <- wa(SumSST ~ ., data = ImbrieKipp, tol.dw = TRUE,
            min.tol = 2, small.tol = "min")
ik.wa

## load V12.122 core data
data(V12.122)
V12.122 <- V12.122 / 100

## predict summer sea-surface temperature for V12.122 core
set.seed(2)
v12.pred <- predict(ik.wa, V12.122, CV = "bootstrap", n.boot = 100)

## draw the fitted reconstruction
reconPlot(v12.pred, use.labels = TRUE, display = "bars")

## extract the model performance stats
performance(v12.pred)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
