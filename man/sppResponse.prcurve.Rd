\name{sppResponse}
\alias{sppResponse}
\alias{sppResponse.prcurve}

\title{Species responses along gradients.}
\description{The fitted responses of species along gradients are
  estimated or extracted from appropriate objects.}
\usage{
sppResponse(x, ...)

\method{sppResponse}{prcurve}(x, n = 100, ...)
}
\arguments{
  \item{x}{an R object.}
  \item{n}{numeric; the number of locations on the gradient to evaluate
    the response curve.}
  \item{\dots}{additional arguments passed to other methods.}
}
\details{
  \code{sppResponse} estimates species responses along indicated
  gradients.

  There is currently no \code{"default"} method and the only specified
  method supplied is for objects fitted by \code{\link{prcurve}}. This
  method extracts the fitted responses of species along the principal
  curve and is a useful diagnostic for identifying overly-complex
  curves.
}
\value{
  A list is returned with components \code{observed} and
  \code{fitted.values} containing the observed and fitted values of the
  species response and gradient respectively. Each is a list with two
  components, \code{gradient} and \code{response}, containing the
  gradient and response values.
}
\author{Gavin L. Simpson}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{prcurve}} for one function that can be used with
  \code{sppResponse}. A \code{\link{plot}} method is available; see
  \code{\link{plot.sppResponse}} for details.
}
\examples{
## Load the Abernethy Forest data set
data(abernethy)

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## Fit the principal curve using varying complexity of smoothers
## for each species
aber.pc <- prcurve(abernethy2, method = "ca", trace = TRUE,
                   vary = TRUE, penalty = 1.4)

## Extract the fitted species response curves
resp <- sppResponse(aber.pc)

## Look at only the most abundant/frequently occurring taxa
take <- chooseTaxa(abernethy2, max.abun = 25, n.occ = 10, value = FALSE)
layout(matrix(1:12, ncol = 3))  	# split device into panels 
plot(resp, which = take)
layout(1)				# reset device
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
