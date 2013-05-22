\name{rankDC}
\alias{rankDC}
\alias{print.rankDC}
\alias{plot.rankDC}
\alias{dotplot.rankDC}

\title{Rank correlation between environmental and species dissimilarities.}
\description{
  Computes the rank correlation between distances between sites in terms
  of gradient variables, such as environmental ones, and species
  composition.
}
\usage{
rankDC(g, y, dc = c("chord", "bray", "euclidean", "chi.square", "gower"),
       method = "spearman")

\method{plot}{rankDC}(x, sort = TRUE, decreasing = FALSE,
     xlab = "Rank correlation", color = "blue",
     pch = 21, bg = "blue", lcolor = "grey",
     lty = "solid", \dots)

\method{dotplot}{rankDC}(x, data = NULL, sort = TRUE, decreasing = FALSE,
     xlab = "Rank correlation", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{the gradient values; usually a data frame of environmental
    data.}
  \item{y}{the species data; usually a data frame.}
  \item{dc}{character; the set of dissimilarity coefficients for which
    rank correlations with gradient distance are to be computed.}
  \item{method}{character; the correlation method to use. See the
    \code{method} argument of \code{\link{cor}}.}
  \item{x}{an object of class \code{"rankDC"}.}
  \item{data}{NULL; ignored, only present for purposes of conformance to
     generic definition.}
  \item{sort, decreasing}{logical; should observations be sorted prior
     to plotting? If so, should they be sorted in order of decreasing
     size?}
  \item{xlab}{The x-axis label for the plot.}
  \item{color, pch, bg, lcolor, lty}{arguments passed to
     \code{\link{dotchart}}.}
  \item{\dots}{arguments passed to other methods, including
     \code{\link{dotchart}} and \code{\link[lattice]{dotplot}}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  A named vector of rank correlations is returned.
}
\author{Gavin L. Simpson, based on \code{\link{rankindex}} from the
  vegan package.}

\seealso{
  \code{\link{rankindex}} from package vegan.
}
\examples{
data(swappH)
data(swapdiat)

rc <- rankDC(swappH, swapdiat, dc = c("chord","euclidean","gower"))

## base plot uses dotchart()
plot(rc)

## lattice version of the base plot
dotplot(rc)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ methods }
\keyword{ utilities }
\keyword{ multivariate }
