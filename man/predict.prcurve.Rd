\name{predict.prcurve}
\alias{predict.prcurve}
\alias{fitted.prcurve}

\title{Predict new locations & fitted values on a principal curve}
\description{
  Locations on a fitted principal curve are predicted by projecting the
  new observations in \eqn{m} dimensions on to the corresponding closest
  point on the curve. Fitted values for data used to fit the curve are
  available with respect to the principal curve or to the individual
  smooth functions.
}
\usage{
\method{predict}{prcurve}(object, newdata, \dots)

\method{fitted}{prcurve}(object, type = c("curve","smooths"), \dots)
}

\arguments{
  \item{object}{
    an object of class \code{\link{prcurve}}.
  }
  \item{newdata}{
    a matrix or data frame of new observations within the space of the
    orginal data. Variables are matched against those of the original
    data via their \code{names} or \code{colnames}. If a data frame is
    supplied, it is converted to a matrix via \code{\link{data.matrix}}.
  }
  \item{type}{
    character; the type of fitted values to return.
  }
  \item{\dots}{
    other arguments passed to other methods. Not currently used.
  }
}
\details{
  Fitting a principal curve involves two procedures. In one, the current
  curve is bent towards the data via the fitting of spline functions
  with distance along the curve as the predictor variable and each
  variable in turn as the response. The second procedure, a projection
  step, involves projecting the observed points in \eqn{m} dimensions on
  to locations along the current curve to which they are closest in the
  hyperspace.

  Given a fitted curve, the projection step can be used to find new
  points on the fitted curve by projecting the new points located in the
  hyperspace on to points on the curve to which they are closest.

  Fitted values are available for the data used to the fit the principal
  curve. There are two types of fitted value available. For \code{type =
  "curve"}, the fitted locations on the principal curve. For \code{type
  = "smooths"}, the fitted values of the variables from the individual
  smooth functions with respect to distance along the principal curve.
}
\value{
  A matrix of points in the space of the original data. Rows correspond
  to the new samples and columns to the variables (ordered as per the
  original data used to fit the curve).

  How these points are ordered along the fitted curve is contained in
  attributed \code{tag}.
}
\author{
  Gavin L. Simpson
}
\seealso{
  See \code{\link{prcurve}} for details on fitting principal curves and
  an example.
}

\keyword{ methods }
