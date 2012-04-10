\name{wa}
\alias{wa}
\alias{wa.default}
\alias{wa.formula}
\alias{print.wa}
\alias{fitted.wa}
\alias{residuals.wa}
\alias{coef.wa}
\alias{waFit}
\title{Weighted averaging transfer functions}
\description{
  Implements the weighted averaging transfer function
  methodology. Tolerance down-weighting and inverse and classicial
  deshrinking are supported.
}
\usage{
wa(x, \dots)

\method{wa}{default}(x, env,
   deshrink = c("inverse", "classical", "expanded", "none"),
   tol.dw = FALSE, useN2 = TRUE,
   na.tol = c("min","mean","max"),
   small.tol = c("min","fraction","absolute"),
   min.tol = NULL, f = 0.1, ...)

\method{wa}{formula}(formula, data, subset, na.action,
   deshrink = c("inverse", "classical", "expanded", "none"),
   tol.dw = FALSE, useN2 = TRUE, na.tol = c("min","mean","max"),
   small.tol = c("min","fraction","absolute"), min.tol = NULL,
   f = 0.1,..., model = FALSE)

\method{fitted}{wa}(object, \dots)

\method{residuals}{wa}(object, \dots)

\method{coef}{wa}(object, \dots)

waFit(x, y, tol.dw, useN2, deshrink, na.tol, small.tol,
      min.tol, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The species training set data}
  \item{env, y}{The response vector}
  \item{deshrink}{Which deshrinking method to use? One of
    \code{"inverse"} or \code{"classical"}, \code{"expanded"} or
    \code{"none"}}
  \item{tol.dw}{logical; should species with wider tolerances be given
    lower weight?}
  \item{useN2}{logical; should Hill's N2 values be used to produce
    un-biased tolerances?}
  \item{na.tol}{character; method to use to replace missing (\code{NA})
    tolerances in WA computations. Missing values are replaced with the
    minimum, average or maximum tolerance observed that is not
    missing.}
  \item{small.tol}{character; method to replace small tolerances. See Details.}
  \item{min.tol}{numeric; threshold below which tolerances are treated
    as being \sQuote{small}. Default is not to replace small tolerances.}
  \item{f}{numeric, \eqn{0 < f < 1}; fraction of environmental gradient
    \code{env} to replace small tolerances with if \code{small.tol =
      "fraction"} is specified.}
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
  where \code{response} is the (numeric) response vector (the variable to
  be predicted) and \code{terms} is a series of terms which specifies a
  linear predictor for \code{response}.  A terms specification of the
  form \code{first + second} indicates all the terms in \code{first}
  together with all the terms in \code{second} with duplicates
  removed. A specification of \code{.} is shorthand for all terms in
  \code{data} not already included in the model.

  Species that have very small tolerances can dominate reconstructed
  values if tolerance down-weighting is used. In \code{wa}, small
  tolerances are defined as a tolerance that is \eqn{<}
  \code{min.tol}. The default is to not replace small tolerances, and
  the user needs to specify suitable values of \code{min.tol}. Function
  \code{\link{tolerance}} may be of use in computing tolerances before
  fitting the WA model.

  Small tolerances can be adjusted in several ways:

  \describe{
    \item{\code{min} }{small tolerances are replaced by the smallest
      observed tolerance that is greater than, or equal to,
      \code{min.tol}. With this method, the replaced values will be no
      smaller than any other observed tolerance. This is the default in
      \pkg{analogue}.}
    \item{\code{fraction} }{small tolerances are replaced by the
      fraction, \code{f}, of the observed environmental gradient in the
      training set, \code{env}.}
    \item{\code{absolute} }{small tolerances are replaced by
      \code{min.tol}.}
  }

  Function \code{waFit} is the workhorse implementing the actual WA
  computations. It performs no checks on the input data and returns a
  simple list containing the optima, tolernances, model tolerances,
  fitted values, coefficients and the numbers of samples and
  species. See Value below for details of each component.
}
\value{
  An object of class \code{"wa"}, a list with the following components:
  
  \item{wa.optima}{The WA optima for each species in the model.}
  \item{tolerances}{The actual tolerances calculated (these are weighted
    standard deviations).}
  \item{model.tol}{The tolerances used in the WA model
    computations. These will be similar to \code{tol}, but will no
    contain any \code{NA}s and any small tolerances will have been
    replaced with the appropriate value.}
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
  \item{options.tol}{A list, containing the values of the arguments
  \code{useN2}, \code{na.tol}, \code{small.tol}, \code{min.tol}, and
  \code{f}.}
  \item{terms, model}{Model \code{\link{terms}} and \code{\link{model.frame}}
    components. Only returned by the \code{formula} method of
    \code{wa}.}
}
%\references{ }
\author{Gavin L. Simpson and Jari Oksanen}
\seealso{\code{\link{mat}} for an alternative transfer function method.}
\examples{
data(ImbrieKipp)
data(SumSST)

## fit the WA model
mod <- wa(SumSST ~., data = ImbrieKipp)
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

## caterpillar plot of optima and tolerances
caterpillarPlot(mod)                 ## observed tolerances
caterpillarPlot(mod, type = "model") ## with tolerances used in WA model

## tolerance DW
mod2 <- wa(SumSST ~ ., data = ImbrieKipp, tol.dw = TRUE,
           min.tol = 2, small.tol = "min")
mod2

## compare actual tolerances to working values
with(mod2, rbind(tolerances, model.tol))
}
\keyword{methods}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
