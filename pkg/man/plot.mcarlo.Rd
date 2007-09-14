\name{plot.mcarlo}
\alias{plot.mcarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting of Monte Carlo simulated dissimilarity distributions}
\description{
  A \code{\link{plot.lm}}-like plotting function for objects of class
  \code{"mcarlo"} to visualise the simulated distribution of
  dissimilarities.
}
\usage{
\method{plot}{mcarlo}(x,
     which = c(1:2),
     alpha = 0.05,
     caption = c("Distribution of dissimilarities",
       expression(paste("Simulated probability Pr (Dissim < ",
           alpha, ")"))),
     col.poly = "lightgrey",
     border.poly = "lightgrey",
     ask = prod(par("mfcol")) < length(which) &&
                                  dev.interactive(),
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"mcarlo"}, usually the result of a
    call to \code{\link{mcarlo}}.}
  \item{which}{numeric; which of the plots should be produced?}
  \item{alpha}{numeric; the Monte Carlo significance level to be marked
    on the cumulative frequency plot.}
  \item{caption}{character, length 2; captions to appear above the
    plots.}
  \item{col.poly, border.poly}{character; the colour to draw the region
    and border of the polygon enclosing the Monte Carlo significance on
    the cummulative frequency plot.} 
  \item{ask}{logical; should the function wait for user confirmation to
    draw each plot? If missing, the function makes a reasonable attempt
    to guess the current situation and act accordingly.}
  \item{\dots}{additional graphical parameters to be passed to the
    plotting functions. Currently ignored.}
}
\details{
  The "Distribution of dissimilarities" plot produces a histogram and
  kernel density estimate of the distribution of simulated dissimilarity
  values.

  The "Simulated probability" plot shows a cumulative probability
  function of the simulated dissimlarity values, and highlights the
  proportion of the curve that is less than \code{alpha}.
}
\value{
  One or more plots on the current device.
}
\references{

  Sawada, M., Viau, A.E., Vettoretti, G., Peltier, W.R. and Gajewski,
  K. (2004) Comparison of North-American pollen-based temperature and
  global lake-status with CCCma AGCM2 output at 6 ka. \emph{Quaternary
    Science Reviews} \strong{23}, 87--108.

}
\author{Gavin L. Simpson}
\seealso{\code{\link{mcarlo}}}
\examples{
## continue the RLGH example from ?join
example(join)

## perform the modified method of Sawada (2004) - paired sampling,
## with replacement
swap.mcarlo <- mcarlo(swapdiat, method = "SQchord", nsamp = 1000,
                      type = "paired", replace = FALSE)
swap.mcarlo

## plot the simulated distribution
par(mfrow = c(2,1))
plot(swap.mcarlo)
par(mfrow = c(1,1))
}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
