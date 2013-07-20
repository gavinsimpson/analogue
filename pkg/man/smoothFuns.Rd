\name{smoothers}
\alias{smoothSpline}
\alias{smoothGAM}

\title{
  Smoother plugin function for use in fitting a principal curve
}

\description{
  Functions to be used as plugins to \code{\link{prcurve}} that fit
  smooth functions to each variable that, when combined, give the
  principal curve. The functions act as wrappers to the main fitting
  functions, which currently include \code{\link{smooth.spline}} and
  \code{\link[mgcv]{gam}}.
}

\usage{
smoothSpline(lambda, x, choose = TRUE, complexity, ...,
             penalty = 1, cv = FALSE, keep.data = FALSE,
             control.spar = list(low = 0))

smoothGAM(lambda, x, choose = TRUE, complexity, bs = "tp", ...,
          family = gaussian(), method = "REML", select = FALSE,
          control = list())
}
\arguments{
  \item{lambda}{the current projection function; the position that each
    sample projects to on the current principal curve. This is the
    predictor variable or covariate in the smooth function.}
  \item{x}{numeric vector; used as the response variable in the smooth
    function. The principal curve algorithm fits a separate scatterplot
    smoother (or similar smoother) to each  variable in \code{X}
    in turn as the response.}
  \item{choose}{logical; should the underlying smoother function be
    allowed to choose the degree of smooth complexity for each
    variable?}
  \item{complexity}{numeric; the complexity of the fitted smooth
    functions.}
  \item{penalty, cv, keep.data, control.spar}{arguments to
    \code{\link{smooth.spline}}.}
  \item{bs, family}{arguments to \code{\link[mgcv]{s}}.}
  \item{method, select, control}{arguments to \code{\link{gam}}.}
  \item{\dots}{arguments passed on the the underlying function
    \code{\link{smooth.spline}} and users should read that function's
    help page for further details.}
}

\value{
  An object of class \code{"prcurveSmoother"} with the following
  components:

  \item{lambda}{for each observations, its arc-length from the beginning
    of the curve.}
  \item{x}{numeric vector of response values.}
  \item{fitted.values}{numeric vector of fitted values for the
    observations generated from the fitted smooth function.}
  \item{complexity}{numeric; the degrees of freedom used for the smooth
    function. The exact details of what these pertain to are in the help
    for the respective fitting functions \code{\link{smooth.spline}} and
    \code{\link[mgcv]{gam}}.}
  \item{model}{the object fitted by the wrapped fitting function.}
}

\author{
  Gavin L. Simpson
}

\seealso{
  \code{\link{prcurve}} for how these functions are used.
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{smooth}