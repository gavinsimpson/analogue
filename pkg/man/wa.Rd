\name{wa}
\alias{wa}
\alias{wa.default}
\alias{wa.formula}
\alias{print.wa}
\alias{fitted.wa}
\alias{residuals.wa}
\alias{coef.wa}
\title{Weighted averaging transfer functions}
\description{
Implements the weighted averaging transfer function
methodology. Tolerance down-weighting and inverse and classicial
deshrinking are supported.
}
\usage{
wa(x, \dots)

\method{wa}{default}(x, env, deshrink = c("inverse", "classical", "expanded", "none"),
   tol.dw = FALSE, ...)

\method{wa}{formula}(formula, data, subset, na.action,
   deshrink = c("inverse", "classical", "expanded", "none"),
   tol.dw = FALSE, model = FALSE, ...)

\method{fitted}{wa}(object, \dots)

\method{residuals}{wa}(object, \dots)

\method{coef}{wa}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The species training set data}
  \item{env}{The response vector}
  \item{deshrink}{Which deshrinking method to use? One of
    \code{"inverse"} or \code{"classical"}, \code{"expanded"} or
    \code{"none"}}
  \item{tol.dw}{logical; should species with wider tolerances be given
    lower weight? Currently ignored as tolerance down weighting is not
    yet implemented.}
  \item{formula}{a model formula}
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables specified on the RHS of the model formula. If not found in
    \code{data}, the  variables are taken from
    \code{environment(formula)}, typically the environment from which
    \code{wa} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  The default is set by the
    \code{na.action} setting of \code{options}, and is \code{na.fail} if
    that is unset.  The 'factory-fresh' default is \code{na.omit}.
    Another possible value is \code{NULL}, no action. Value
    \code{na.exclude} can be useful.}
  \item{model}{logical. If \code{TRUE} the model frame is returned.}
  \item{object}{an Object of class \code{"wa"}, the result of a call to
    \code{wa}.}
  \item{\dots}{arguments to other methods.}
}
\details{
  A typical model has the form \code{response ~ terms}
  where\code{response} is the (numeric) response vector (the variable to
  be predicted) and \code{terms} is a series of terms which specifies a
  linear predictor for \code{response}.  A terms specification of the
  form \code{first + second} indicates all the terms in \code{first}
  together with all the terms in \code{second} with duplicates
  removed. A specification of \code{.} is shorthand for all terms in
  \code{data} not already included in the model.
}
\value{
  An object of class \code{"wa"}, a list with the following components:
  
  \item{wa.optima}{The WA optima for each species in the model.}
  \item{fitted.values}{The fitted values of the response for each of the
    training set samples.}
  \item{residuals}{Model residuals.}
  \item{coefficients}{Deshrinking coefficients.}
  \item{rmse}{The RMSE of the model.}
  \item{r.squared}{The coefficient of determination of the observed and
    fitted values of the response.}
  \item{avg.bias, max.bias}{The average and maximum bias statistics.}
  \item{n.samp, n.spp}{The number of samples and species in the training
    set.}
  \item{deshrink}{The deshrinking regression method used.}
  \item{tol.dw}{logical; was tolerance down-weighting applied?}
  \item{call}{The matched function call.}
  \item{orig.x}{The training set species data.}
  \item{orig.env}{The response data for the training set.}
  \item{terms, model}{Model \code{\link{terms}} and \code{\link{model.frame}}
    components. Only returned by the \code{formula} method of
    \code{wa}.}
}
%\references{ }
\author{Gavin L. Simpson and Jari Oksanen}
\seealso{\code{\link{mat}} for an alternative transfer function method.}
\examples{
data(swapdiat)
data(swappH)
swapdiat <- swapdiat / 100

## fit the WA model
mod <- wa(swappH ~., data = swapdiat)
mod

## extract the fitted values
fitted(mod)

## residuals for the training set
residuals(mod)

## deshrinking coefficients
coef(mod)

## diagnostics plots
par(mfrow = c(1,2))
plot(mod)
par(mfrow = c(1,1))
}
\keyword{methods}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
