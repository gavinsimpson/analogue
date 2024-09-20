\name{plot.prcurve}
\alias{plot.prcurve}
\alias{lines.prcurve}
\title{
Plot a fitted principal curve in PCA space
}
\description{
Projects the principal curve into PCA space and draws it and the
underlying data in a biplot.
}
\usage{
\method{plot}{prcurve}(x, axes = 1:2, scaling = 0, segments = TRUE,
     col = "red", col.seg = "forestgreen", lwd = 2,
     lwd.seg = 1, ...)

\method{lines}{prcurve}(x, axes = 1:2, scaling = 0, segments = TRUE,
      col = "red", col.seg = "forestgreen", lwd = 2,
      lwd.seg = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"prcurve"}.}
  \item{axes}{numeric vector of length 2; this is passed to the
    \code{choices} argument of the \code{\link[vegan]{scores}}
    function.}
  \item{scaling}{numeric; the scaling to use. See
    \code{\link[vegan]{scores.rda}} for the available options. The default is
    not to scale the scores, but \code{scaling = 1} might be a useful
    alternative.}
  \item{segments}{logical; should segments be drawn between the observed
    points to the location on the principal curve on to which they
    project.}
  \item{col}{The colour to draw the principal curve in.}
  \item{col.seg}{The colour to draw the segments in.}
  \item{lwd, lwd.seg}{The line thickness used to draw the
      principal curve and segments respectively.}
  \item{\dots}{additional arguments passed on to \code{points} when
    drawing the observations in PCA space.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A plot on the currently active device. The function does not return
anything.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Gavin L. Simpson
}
%% \note{
%%   Note that all plotting is done using \code{scaling == 0} at the
%%   moment, and as such the \code{lines()} method only makes sense when
%%   added to an underlying PCA in the same scaling. See the Examples.
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{prcurve}}; \code{\link[vegan]{rda}} for the code used to perform
the PCA.
}
\examples{
## Load the Abernethy Forest data
data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit the principal curve using varying complexity of smoothers
## for each species
aber.pc <- prcurve(abernethy2, method = "ca", trace = TRUE,
                   vary = TRUE, penalty = 1.4)

## Plot the curve
plot(aber.pc)

## The lines() method can be used to add the principal curve to an
## existing plot
ord <- rda(abernethy2)
plot(ord, scaling = 1)
lines(aber.pc, scaling = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
