\name{plot.prcurve}
\alias{plot.prcurve}
\title{
Plot a fitted principal curve in PCA space
}
\description{
Projects the principal curve into PCA space and draws it and the
underlying data in a biplot.
}
\usage{
\method{plot}{prcurve}(x, data, axes = 1:2, seg = TRUE, col.seg = "forestgreen",
     col.curve = "red", lwd.curve = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"prcurve"}.}
  \item{data}{The data the principal curve was fitted to.}
  \item{axes}{numeric vector of length 2; this is passed to the
    \code{choices} argument of the \code{\link[vegan]{scores}} function.}
  \item{seg}{logical; should segments be drawn between the observed
    points to the location on the principal curve on to which they
    project.} 
  \item{col.seg}{The colour to draw the segments in.}
  \item{col.curve}{The colour to draw the principal curve in.}
  \item{lwd.curve}{The line thickness used to draw the principal curve.}
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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{prcurve}}; \code{\link{rda}} for the code used to perform
the PCA.
}
\examples{
## Load the Abernethy Forest data
data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit the principal curve using varying complexity of smoothers
## for each species
aber.pc2 <- prcurve(abernethy2, method = "ca", trace = TRUE,
                    vary = TRUE, penalty = 1.4)

## Plot the curve
plot(aber.pc2, abernethy2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
