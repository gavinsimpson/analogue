\name{bootstrap}
\alias{bootstrap}
\alias{bootstrap.default}
\alias{bootstrap.mat}
\alias{print.bootstrap.mat}
\alias{residuals.bootstrap.mat}
\alias{resid.bootstrap.mat}
\alias{print.residuals.bootstrap.mat}
\alias{fitted.bootstrap.mat}
\alias{print.fitted.bootstrap.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap estimation and errors }
\description{
  Function to calculate bootstrap statistics for transfer function
  models such as bootstrap estimates, model RMSEP, sample specific
  errors for predictions and summary statistics such as bias and
  \eqn{R^2} between oberved and estimated
  environment.

  \code{\link{residuals}} method for objects of class
  \code{"bootstrap.mat"}.
}
\usage{

bootstrap(object, \dots)

\method{bootstrap}{default}(object, \dots)

\method{bootstrap}{mat}(object, newdata, newenv, k,
          weighted = FALSE, n.boot = 1000, \dots)

\method{fitted}{bootstrap.mat}(object, k, \dots)

\method{residuals}{bootstrap.mat}(object, which = c("model", "bootstrap"), \dots) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object of class \code{"mat"} for which bootstrap
    statistics are to be generated, or an object of class
    \code{"bootstrap.mat"} from which fitted values or residuals are
    extracted.}
  \item{newdata}{a data frame containing samples for which bootstrap
    predictions and sample specific errors are to be generated. May be
    missing --- See Details. \code{"newdata"} must have the same number
    of columns as the training set data.}
  \item{newenv}{a vector containing environmental data for samples
    in \code{"newdata"}. Used to calculate full suite of errors for new
    data such as a test set with known environmental values. May be
    missing --- See Details. \code{"newenv"} must have the same number
    of rows as \code{"newdata"}.}
  \item{k}{numeric; how many modern analogues to use to generate the
    bootstrap statistics (and, if requested, the predictions), fitted
    values or residuals.}
  \item{weighted}{logical; should the weighted mean of the environment
    for the \code{"k"} modern analogues be used instead of the mean?}
  \item{n.boot}{Number of bootstrap samples to take.}
  \item{which}{character; which set of residuals to return, the
    model residuals or the residuals of the bootstrap-derived
    estimates?}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  \code{bootstrap} is a fairly flexible function, and can be called with
  or without arguments \code{newdata} and \code{newenv}.

  If called with only \code{object} specified, then bootstrap estimates
  for the training set data are returned. In this case, the returned
  object will not include component \code{predictions}.

  If called with both \code{object} and \code{newdata}, then in addition
  to the above, bootstrap estimates for the new samples are also
  calculated and returned. In this case, component \code{predictions}
  will contain the apparent and bootstrap derived predictions and
  sample-specific errors for the new samples.

  If called with \code{object}, \code{newdata} and \code{newenv}, then
  the full \code{bootstrap} object is returned (as described in the
  Value section below). With environmental data now available for the
  new samples, residuals, RMSE(P) and \eqn{R^2} and bias statistics can
  be calculated.

  The individual components of \code{predictions} are the same as those
  described in the components relating to the training set data. For
  example, \code{returned.object$predictions$bootstrap} contains the
  components as \code{returned.object$bootstrap}.

  It is not usual for environmental data to be available for the new
  samples for which predictions are required. In normal
  palaeolimnological studies, it is more likely that \code{newenv} will
  not be available as we are dealing with sediment core samples from the
  past for which environmental data are not available. However, if
  sufficient training set samples are available to justify producing a
  training and a test set, then \code{newenv} will be available, and
  \code{bootstrap} can accomodate this extra information and calculate
  apparent and bootstrap estimates for the test set, allowing an
  independent assessment of the RMSEP of the model to be performed.
  
  Typical usage of \code{residuals} is
  \preformatted{
    resid(object, which = c("model", "bootstrap"), \dots)}
}
\value{

  For \code{bootstrap.mat} an object of class \code{"bootstrap.mat"} is
  returned. This is a complex object with many components and is
  described in \code{\link{bootstrapObject}}.
  
  For \code{residuals}, a list containg the requested residuals and
  metadata, with the following components:
  \item{model}{Leave one out residuals for the MAT-estimated model.}
  \item{bootstrap}{residuals for the bootstrapped MAT model.}
  \item{k}{numeric; indicating the size of model used in estimates and
    predictions.}
  \item{n.boot}{numeric; the number of bootstrap samples taken.}
  \item{auto}{logical; whether \code{"k"} was choosen automatically or
    user-selected.}
  \item{weighted}{logical; whether the weighted mean was used instead of
    the mean of the environment for \emph{k}-closest analogues.}
}
\references{

  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C. and ter Braak,
  C.J.F. (1990). Diatoms and pH reconstruction. \emph{Philosophical
    Transactions of the Royal Society of London; Series B}, \strong{327};
  263--278.

}
\author{Gavin L. Simpson }
\seealso{\code{mat}, \code{plot.mat}, \code{summary.bootstrap.mat},
  \code{\link{residuals}}}
\examples{
\dontshow{op <- options(digits = 3)}
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

## Imbrie and Kipp foraminfera sea-surface temperature 
## fit the MAT model using the squared chord distance measure
ik.mat <- mat(ImbrieKipp, SumSST, method = "SQchord")

## bootstrap training set
ik.boot <- bootstrap(ik.mat, n.boot = 100)
ik.boot
summary(ik.boot)

## Bootstrap fitted values for training set
fitted(ik.boot)

## residuals
resid(ik.boot) # uses abbreviated form
\dontshow{options(op)}
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{methods}
