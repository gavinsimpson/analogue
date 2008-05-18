\name{tran}
\alias{tran}
\concept{transformation}
\title{Common data transformations and standardizations}
\description{
  Provides common data transformations and standardizations useful for
  palaeoecological data. The function acts as a wrapper to function
  \code{\link[vegan]{decostand}} in package vegan for several of the
  available options.
}
\usage{
tran(x, method, a = 1, b = 0, base = exp(1), na.rm = FALSE,
     na.value = 0, ...)
}
\arguments{
  \item{x}{A matrix-like object.}
  \item{method}{transformation or standardization method to apply. See
    Details for available options.}
  \item{a}{Constant to multiply \code{x} by. \code{method = "log"}
    only. Can be a vector, in which case the vector of values to
    multiply each column of \code{x} by.}
  \item{b}{Constant to add to \code{x} before taking logs. \code{method
      = "log"} only. Can be a vector, in which case the vector of values
    to add to each column of \code{x}.}
  \item{base}{the base with respect to which logarithms are
    computed. See \code{\link{log}} for further details. The default is
    to compute natural logarithms.}
  \item{na.rm}{Should missing values be removed before some computations?}
  \item{na.value}{The value with which to replace missing values
    (\code{NA}).}
  \item{\dots}{Arguments passed to \code{\link[vegan]{decostand}}.}
}
\details{
  The function offers following standardization methods for community
  data:
  \itemize{
    \item \code{sqrt}: take the square roots of the observed values.
    \item \code{cubert}: take the cube root of the observed values.
    \item \code{log}: take the logarithms of the observed values. The
    tansformation applied can be modified by constants \code{a} and
    \code{b} and the \code{base} of the logarithms. The transformation
    applied is \eqn{x^* = \log_{\mathrm{base}}(ax + b)}{x^* =
      log[base](ax + b).}
    \item \code{reciprocal}: returns the multiplicative inverse or
    reciprocal, \eqn{1/x}, of the observed values.
    \item \code{freq}: divide by column (variable, species)  maximum and
    multiply by the number of non-zero items, so that the average of
    non-zero entries is 1 (Oksanen 1983).
    \item \code{center}: centre all variables to zero mean.
    \item \code{range}: standardize values into range 0 \dots 1. If all
    values are constant, they will be transformed to 0.
    \item \code{percent}: convert observed count values to percentages.
    \item \code{proportion}: convert observed count values to proportions.
    \item \code{standardize}: scale \code{x} to zero mean and unit
    variance.
    \item \code{pa}: scale \code{x} to presence/absence scale (0/1).
    \item \code{missing}: replace missing values with \code{na.value}.
    \item \code{chi.square}: divide by row sums and square root of
    column sums, and adjust for square root of matrix total
    (Legendre & Gallagher 2001). When used with the Euclidean
    distance, the distances should be similar to the the
    Chi-square distance used in correspondence analysis. However, the
    results from \code{\link{cmdscale}} would still differ, since
    CA is a weighted ordination method.
    \item \code{hellinger}: square root of observed values that have
    first been divided by row (site) sums (Legendre & Gallagher 2001).
    \item \code{wisconsin}: applies the Wisconsin double
    standardization, where columns (species, variables) are first
    standardized by maxima and then sites (rows) by site totals.
    \item \code{pcent2prop}: convert percentages to proportions.
    \item \code{prop2pcent}: convert proportions to percentages.
  }
}
\value{
  Returns the suitably transformed or standardized \code{x}. If \code{x}
  is a data frame, the returned value is like-wise a data frame. The
  returned object also has an attribute \code{"tran"} giving the name of
  applied transformation or standardization \code{"method"}.
}
\references{
  Legendre, P. & Gallagher, E.D. (2001) Ecologically meaningful
  transformations for ordination of species data. \emph{Oecologia}
  \strong{129}; 271-280.
  
  Oksanen, J. (1983) Ordination of boreal heath-like vegetation with
  principal component analysis, correspondence analysis and
  multidimensional scaling. \emph{Vegetatio} \strong{52}; 181-189.
}
\author{Gavin L. Simpson. Much of the functionality of \code{tran} is
  provided by \code{\link[vegan]{decostand}}, written by Jari Oksanen.}
\seealso{\code{\link[vegan]{decostand}}}
\examples{
data(swapdiat)
## convert percentages to proportions
sptrans <- tran(swapdiat, "pcent2prop")

## apply Hellinger transformation
spHell <- tran(swapdiat, "hellinger")
}
\keyword{ multivariate }
\keyword{ manip }
