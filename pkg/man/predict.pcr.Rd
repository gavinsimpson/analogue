\name{predict.pcr}
\alias{predict.pcr}

\title{Predicted values from a principal components regression}

\description{
  Calculates predicted values from a fitted principal components
  regression model. Leave-one-out, bootstrap of n-fold crossvalidated
  predictions are also intended (but not yet implemented).
}

\usage{
\method{predict}{pcr}(object, newdata, ncomp = seq_along(object$ncomp),
        CV = c("none", "LOO", "bootstrap", "nfold"),
        verbose = FALSE, nboot = 100, nfold = 5, \dots)
}

\arguments{
  \item{object}{a fitted model of class \code{"pcr"}, the result of a
    call to \code{\link{pcr}}.}
  \item{newdata}{data frame of new observations for which predictions
    are sought.}
  \item{ncomp}{numeric; the PCR components for which predictions are
    sought. Can be a vector in which case predictions for multiple
    components are computed.}
  \item{CV}{character; the type of crossvalidation required. Currently,
    no crossvalidation methods are implemented.}
  \item{verbose}{logical; should progress on crossvalidation be printed
    to the console?}
  \item{nboot}{numeric; the number of bootstrap samples to draw, or in
    the case of \code{CV = "nfold"} the number of repeats of n-fold CV
    to perform.}
  \item{nfold}{numeric; the number of folds to split data into.}
  \item{\dots}{arguments passed to other methods.}
}

\details{
  \code{predict.pcr} arranges for any transformation applied to the
  training data to be applied to the \code{newdata} prior to
  prediction.
}

\value{
  A matrix of predicted values with rows representing samples in
  \code{newdata} and columns, the PCR components requested via
  \code{ncomp}.
}

\author{Gavin L. Simpson}

\seealso{\code{\link{pcr}}}

\examples{
## Load the Imbrie & Kipp data and
## summer sea-surface temperatures
data(ImbrieKipp)
data(SumSST)

## choose 10 samples to act as a test set, for illustration
take <- c(5,58,31,51,42,28,30,57,8,50)

## normal interface and apply Hellinger transformation
mod <- pcr(ImbrieKipp[-take, ], SumSST[-take], tranFun = Hellinger)

## predictions
predict(mod, ImbrieKipp[take, ], ncomp = 1:4)

}
\keyword{methods}