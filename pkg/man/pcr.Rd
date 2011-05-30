\name{pcr}
\alias{pcr}
\alias{pcr.default}
\alias{pcr.formula}
\alias{print.pcr}
\alias{Hellinger}
\alias{ChiSquare}
\title{Prinicpal component regression transfer function models}
\description{
  Fits a palaeoecological transfer function model using principal
  component regression, using an optional transformation of the matrix
  of predictor variables when these are species abundance data.
}
\usage{

\method{pcr}{default}(x, y, ncomp, tranFun, ...)

\method{pcr}{formula}(formula, data, subset, na.action, ..., model = FALSE)

Hellinger(x, apply = FALSE)

ChiSquare(x, apply = FALSE)
}
\arguments{
  \item{x}{Matrix or data frame of predictor variables. Usually species
    composition or abundance data for transfer function models.}
  \item{y}{Numeric vector; the response variable to be modelled.}
  \item{ncomp}{numeric; number of principal components to build models
    for. If not supplied the largest possible number of components is
    determined.}
  \item{tranFun}{function; a function or name of a function that
    performs a transformation of the predictor variables \code{x}. The
    function must be self-contained as no arguments are passed to the
    function when it is applied. See Details for more information.}
  \item{formula}{a model formula.}
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables specified on the RHS of the model formula. If not found in
    \code{data}, the  variables are taken from
    \code{environment(formula)}, typically the environment from which
    \code{pcr} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  The default is set by the
    \code{na.action} setting of \code{options}, and is \code{na.fail} if
    that is unset.  The 'factory-fresh' default is \code{na.omit}.
    Another possible value is \code{NULL}, no action. Value
    \code{na.exclude} can be useful.}
  \item{model}{logical; If \code{TRUE} the model frame is returned?}
  \item{apply}{logical; should an existing tranformation, using
    pre-computed meta-parameters, be applied?}
  \item{\dots}{Arguments passed to other methods.}
}
\details{
  When applying cross-validation (CV) to transfer function models, any
  transformation of the predictors must be applied separately during
  each iteration of the CV procedure to the part of the data used in
  fitting the model. In the same way, any samples to be predicted from
  the model must use any meta-parameters derived from the training data
  only. For examle, centring is appled to the training data only and the
  variables means used to centre the training data are used to centre
  the test samples. The variable means should not be computed on a
  combination of the training and test samples.

  When using PCR, we might wish to apply a transformation to the species
  data predictor variables such that the PCA of those data preserves a
  dissimilarity coefficient other than the Euclidean distance. This
  transformation is applied to allow PCA to better describe patterns in
  the species data (Legendre & Gallagher 2001).

  How this is handled in \code{pcr} is to take a user-supplied function
  that takes a single argument, the matrix of predictor variables. The
  function should return a matrix of the same dimension as the input. If
  any meta-parameters are required for subsequent use in prediction,
  these should be returned as attribute \code{"parms"}, attached to the
  matrix.

  Two example transformation functions are provided implementing the
  Hellinger and Chi Square transformations of Legendre & Gallagher
  (2001). Users can base their transformation functions on
  these. \code{ChiSquare()} illustrates how meta-parameters should be
  returned as the attribute \code{"parms"}.
}
\value{
  Returns an object of class \code{"pcr"}, a list with the
  following components:

  \item{fitted.values}{matrix; the PCR estimates of the response. The
    columns contain fitted values using C components, where C is the Cth
    column of the matrix.}
  \item{coefficients}{matrix; regression coefficients for the
    PCR. Columns as per \code{fitted} above.}
  \item{residuals}{matrix; residuals, where the Cth column represents a
    PCR model using C components.}
  \item{scores}{}
  \item{loadings}{}
  \item{Yloadings}{}
  \item{xMeans}{numeric; means of the predictor variables in the
    training data.}
  \item{yMean}{numeric; mean of the response variable in the training
    data.}
  \item{varExpl}{numeric; variance explained by the PCR model. These are
    the squares of the singular values.}
  \item{totvar}{numeric; total variance in the training data}
  \item{call}{the matched call.}
  \item{tranFun}{transformation function used. \code{NA} if none
    supplied/used.}
  \item{tranParms}{list; meta parameters used to computed the
    transformed training data.}
  \item{performance}{data frame; cross-validation performance statistics
    for the model.}
  \item{ncomp}{numeric; number of principal components computed}
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

## normal interface and apply Hellinger transformation
mod <- pcr(ImbrieKipp, SumSST, tranFun = Hellinger)
mod

## formula interface, but as above
mod2 <- pcr(SumSST ~ ., data = ImbrieKipp, tranFun = Hellinger)
mod2

}
\keyword{methods}
