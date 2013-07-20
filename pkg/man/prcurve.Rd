\name{prcurve}
\alias{prcurve}
\alias{initCurve}
\alias{print.prcurve}

\title{
Fits a principal curve to m-dimensional data
}
\description{
  A principal curve is a non-parametric generalisation of the principal
  component and is a curve that passes through the \emph{middle} of a
  cloud of data points for a certain definition of `middle'.
}
\usage{
prcurve(X, method = c("ca", "pca", "random", "user"), start = NULL,
        smoother = smoothSpline, complexity, vary = FALSE,
        maxComp, finalCV = FALSE, axis = 1, rank = FALSE,
        stretch = 2, maxit = 10, trace = FALSE, thresh = 0.001,
        plotit = FALSE, \dots)

initCurve(X, method = c("ca", "pca", "random", "user"), rank = FALSE,
          axis = 1, start)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix-like object containing the variables to which the
    principal curve is to be fitted.}
  \item{method}{character; method to use when initialising the principal
    curve. \code{"ca"} fits a correspondence analysis to \code{X} and uses
    the \code{axis}-th axis scores as the initial curve. \code{"pca"} does
    the same but fits a principal components analysis to
    \code{X}. \code{"random"} produces a random ordering as the initial
    curve.}
  \item{start}{numeric vector specifying the initial curve when
    \code{method = "user"}. Must be of length \code{nrow(X)}.}
  \item{smoother}{function; the choice of smoother used to fit the
    principal curve. Currently, the only options are
    \code{smoothSpline}, which is a wrapper to
    \code{\link{smooth.spline}}, and \code{smoothGAM}, which is a
    wrapper to  \code{\link[mgcv]{gam}}.}
  \item{complexity}{numeric; the complexity of the fitted smooth
    functions.
    
    The function passed as argument \code{smoother} should arrange for
    this argument to be passed on to relevant aspect of the underlying
    smoother. In the case of \code{smoothSpline}, complexity is the
    \code{df} argument of \code{\link{smooth.spline}}.}
  \item{vary}{logical; should the complexity of the smoother fitted to
    each variable in \code{X} be allowed to vary (i.e. to allow a more or
    less smooth function for a particular variable. If \code{FALSE} the
    median complexity over all \emph{m} variables is chosen as the fixed
    complexity for all \emph{m} smooths.}
  \item{maxComp}{numeric; the upper limt on the allowed complexity.}
  \item{finalCV}{logial; should a final fit of the smooth function be
    performed using cross validation?}
  \item{axis}{numeric; the ordinaion axis to use as the initial curve.}
  \item{rank}{logical; should rank position on the gradient be used? Not
    yet implemented.}
  \item{stretch}{numeric; a factor by which the curve can be
    extrapolated when points are projected.  Default is 2 (times the
    last segment length).}
  \item{maxit}{numeric; the maximum number of iterations.}
  \item{trace}{logical; print progress on the iterations be printed to
    the console?}
  \item{thresh}{numeric; convergence threshold on shortest distances to
    the curve. The algorithm is considered to have converged when the
    latest iteration produces a total residual distance to the curve
    that is within \code{thresh} of the value obtained during the
    previous iteration.}
  \item{plotit}{logical; should the fitting process be plotted? If
    \code{TRUE}, then the fitted principal curve and observations in
    \code{X} are plotted in principal component space.}
  \item{\dots}{additional arguments are passed solely on to the function
    \code{smoother}.}
}
%\details{
%TODO
%}
\value{
  An object of class \code{"prcurve"} with the following components:

  \item{s}{a matrix corresponding to \code{X}, giving their projections
    onto the curve.}
  \item{tag}{an index, such that \code{s[tag, ]} is smooth.}
  \item{lambda}{for each point, its arc-length from the beginning of the
    curve.}
  \item{dist}{the sum-of-squared distances from the points to their
    projections.}
  \item{converged}{logical; did the algorithm converge?}
  \item{iter}{numeric; the number of iterations performed.}
  \item{totalDist}{numeric; total sum-of-squared distances.}
  \item{complexity}{numeric vector; the complexity of the smoother
    fitted to each variable in \code{X}.}
  \item{call}{the matched call.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Gavin L. Simpson
}
\note{
  The fitting function uses function \code{\link[princurve]{get.lam}} in
    package \pkg{princurve} to find the projection of the data on to the
    fitted curve. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{smoothGAM}} and \code{\link{smoothSpline}} for the
  wrappers fitting smooth functions to each variable.
}
\examples{
## Load Abernethy Forest data set
data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit the principal curve using the median complexity over
## all species
aber.pc <- prcurve(abernethy2, method = "ca", trace = TRUE,
                   vary = FALSE, penalty = 1.4)

## Fit the principal curve using varying complexity of smoothers
## for each species
aber.pc2 <- prcurve(abernethy2, method = "ca", trace = TRUE,
                    vary = TRUE, penalty = 1.4)

\dontrun{
## Fit principal curve using a GAM - currently slow ~10secs
aber.pc3 <- prcurve(abernethy2, method = "ca", trace = TRUE,
                    vary = TRUE, smoother = smoothGAM, bs = "cr")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{nonparametric}
\keyword{smooth}
