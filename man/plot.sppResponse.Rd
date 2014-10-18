\name{plot.sppResponse}
\alias{plot.sppResponse}

\title{Plot species responses along gradients or latent variables}
\description{
  Observed abundances and fitted response curves along
}
\usage{
\method{plot}{sppResponse}(x,
     which = seq_along(x),
     display = c("observed", "fitted"),
     xlab = "Gradient", ylab = "Abundance",
     main = NULL,
     lcol = "red",
     lwd = 2,
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"sppResponse"}.}
  \item{which}{numeric or logical vector; which species (components of
    \code{x}) to plot.}
  \item{display}{character; plot the observed data or the fitted species
    response curve or both?}
  \item{xlab, ylab, main}{titles for the x and y axis label and the main
    title. If \code{main} is \code{NULL}, suitable main titles are taken
    from \code{x}.}
  \item{lwd}{numeric; width of the line used to draw the fitted species
    response function if drawn.}
  \item{lcol}{the colour to use to draw the fitted species response.}
  \item{\dots}{graphical arguments passed \code{plot} and \code{lines}
    used internally.}
}
\details{
  Currently, no attempt is made to split the device into an appropriate
  number of panels. This is up to the user to decide. See the Examples
  section of \code{\link{sppResponse}} for one way to handle this.
}
\value{
  One or more plots, drawn on the current device.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{sppResponse}} for a complete example using fitted
  responses along principal curves.}
\keyword{hplot}
\keyword{methods}
