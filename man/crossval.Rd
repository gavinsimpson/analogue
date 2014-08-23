\name{crossval}
\alias{crossval}
\alias{crossval.wa}
\alias{crossval.pcr}
\alias{print.crossval}
%\alias{predWA}
%\alias{predWAT}
\title{Cross-validation of palaeoecological transfer function models}
\description{
  Performs leave-one-out, \emph{k}-fold, \emph{n} \emph{k}-fold and
  bootstrap cross-validation of palaeoecological transfer function models.
}
\usage{
crossval(obj, ...)

\method{crossval}{wa}(obj, method = c("LOO","kfold","bootstrap"),
         nboot = 100, nfold = 10, folds = 5,
         verbose = getOption("verbose"), ...)

\method{crossval}{pcr}(obj, method = c("LOO","kfold","bootstrap"),
         ncomp, nboot = 100, nfold = 10, folds = 5,
         verbose = getOption("verbose"), ...)

}
\arguments{
  \item{obj}{A fitted transfer function model. Currently, only objects
    of class \code{\link{wa}} and \code{\link{pcr}} are supported.}
  \item{method}{character; type of cross-validation.}
  \item{ncomp}{numeric; number of components to fit, as in models with
    \code{1:ncomp} components.}
  \item{nboot}{numeric; number of bootstrap samples.}
  \item{nfold}{numeric; number of chunks into which the training data
    are split. The \emph{k} in \emph{k}-fold.}
  \item{folds}{numeric; the number of times \emph{k}-fold CV is
    performed.}
  \item{verbose}{logical; should progress of the CV be displayed?}
  \item{\dots}{Arguments passed to other methods.}
}
%\details{
%
%}
\value{
  Returns an object of class \code{"crossval"}, a list with the
  following components:

  \item{fitted.values}{numeric vector; the cross-validated estimates of
    the response.}
  \item{residuals}{numeric vector; residuals computed from the
    cross-validated estimates of the response.}
  \item{performance}{data frame; cross-validation performance statistics
    for the model.}
  \item{CVparams}{list; parameters holding details of the
    cross-validation process.}
  \item{call}{the matched call.}
}
%\references{TO DO}
\author{Gavin L. Simpson}
%\note{
%}
\seealso{\code{\link{wa}}}
\examples{
## Load the Imbrie & Kipp data and
## summer sea-surface temperatures
data(ImbrieKipp)
data(SumSST)
     
## fit the WA model
mod <- wa(SumSST ~., data = ImbrieKipp)
mod

## Leave one out CV
cv.loo <- crossval(mod)
cv.loo

## k-fold CV (k == 10)
cv.kfold <- crossval(mod, method = "kfold", kfold = 10, folds = 1)
cv.kfold

## n k-fold CV (k == 10, n = 10)
cv.nkfold <- crossval(mod, method = "kfold", kfold = 10, folds = 10)
cv.nkfold

## bootstrap with 100 bootstrap samples
cv.boot <- crossval(mod, method = "bootstrap", nboot = 100)
cv.boot

## extract fitted values and residuals
fitted(cv.boot)
resid(cv.boot)

## Principal Components Regression
mpcr <- pcr(SumSST ~., data = ImbrieKipp, ncomp = 10)
crossval(mpcr, method = "kfold", kfold = 10, folds = 2, ncomp = 10)

crossval(mpcr, method = "bootstrap", nboot = 100, ncomp = 10)
}
\keyword{methods}
