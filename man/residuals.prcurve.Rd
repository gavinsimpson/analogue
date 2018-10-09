\name{residuals.prcurve}
\alias{residuals.prcurve}
\alias{resid.prcurve}

\title{
  Residuals of a principal curve fit.
}
\description{
  Returns various representations of the residuals of a principal curve
  fit.
}
\usage{
\method{residuals}{prcurve}(object, which = c("distance", "raw", "smooths", "pca"),
          ...)
}

\arguments{
  \item{object}{an object of class \code{"prcurve"}, the result of a
    call to \code{\link{prcurve}}.}
  \item{which}{character; the type of residuals to return. See Details.}
  \item{\dots}{arguments passed to other methods. See Details.}
}
\details{
  Various types of residual are available for the principal curve. In a
  departure from the usual convention, which residuals are returned is
  controlled via the \code{which} argument. This is to allow users to
  pass a \code{type} argument to the \code{residuals} method for the
  function used to fit the individual smooth functions when \code{which
    = "smooths"}.

  The types of residuals available are

  \describe{
    \item{\code{"distance"}}{the default residual for a principal
      curve. This residual is taken as the Euclidean distance between each
      observations and the point on the principal curve to which it
      projects, in full multivariate space.}
    \item{\code{"raw"}}{raw residuals are the basis for
      \code{"distance"} residuals, and are the difference between the
      observed and fitted values (position on the curve) for each
      observation in terms of each variable in the data set. These
      residuals are in the form of a matrix with number of observation
      \emph{rows} and number of variables \emph{cols}.}
    \item{\code{"smooths"}}{these residuals are the result of calling
      \code{residuals()} on each of the smooth models fitted to the
      individual variables. See below for further details. A matrix of
      the same dimensions as for \code{which = "raw"} is returned.}
    \item{\code{"pca"}}{similar to the raw residuals, but expressed in
      terms of the principal components of the input data. Hence these
      residuals are the difference between each observation's location
      in PCA space and their corresponding location on the curve.}
  }

  For \code{"smooths"} residuals, what is returned is governed by the
  \code{residuals} method available for the smooth model fitted to the
  individual variables. For principal curves fitted using the
  \code{\link{smoothSpline}} plugin, see
  \code{\link{smooth.spline}}. For principal curves fitted via the
  \code{\link{smoothGAM}} plugin, see
  \code{\link[mgcv]{residuals.gam}}.

  \dots can be used to pass additional arguments to these
  \code{residuals} methods. In particular, the \code{type} argument is
  commonly used to choose which type of residual is returned by the
  specific methods.

  In the case of principal curves fitted using the plugin
  \code{\link{smoothSpline}}, residuals for \code{which = "smooths"} are
  only available if the the additional argument \code{keep.data} was
  specified during fitting via \code{\link{prcurve}}. See the examples
  for an illustration of this usage.
}
\value{
  A vector of residual distances (\code{which = "distance"}) or a matrix
  of residuals (for the other options).
}

\author{Gavin L. Simpson}

\seealso{
  \code{\link{prcurve}} for fitting a principal curve.
}

\examples{
  \dontshow{od <- options(digits = 4)}
  ## Load Abernethy Forest data set
  data(abernethy)

  ## Remove the Depth and Age variables
  abernethy2 <- abernethy[, -(37:38)]
  
  ## Fit the principal curve, preserving the data in the smooth.spline
  ## smooth functions fitted via keep.data = TRUE
  aber.pc <- prcurve(abernethy2, method = "ca", keep.data = TRUE)

  ## default "distance" residuals
  res <- resid(aber.pc)
  head(res)

  ## residuals from the underlying smooth models, also illustrates
  ## how to select specific types of residual from the individual
  ## method using argument 'type'
  res <- resid(aber.pc, which = "smooths", type = "deviance")
  dim(res)
  head(res[, 1:5])		# just show a few species residuals
  \dontshow{options(od)}
}

\keyword{ methods }
