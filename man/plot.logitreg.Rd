\name{plot.logitreg}
\alias{plot.logitreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces plots of analogue logistic regression models}
\description{
  Draws the fitted logistic regression function describing the posterior
  probability that two sites are analogues conditional upon the
  dissimilarity between the two samples. Confidence intervals are also
  computed and displayed if requested.
}
\usage{
\method{plot}{logitreg}(x, group = "all", npred = 100,
     conf.int = 0.9, conf.type = c("none", "polygon", "lines"),
     xlab = expression(D[ij]), ylab = "Pr (A+ | d)",
     rug = TRUE, ticksize = 0.02,
     col = "red", ref.col = "lightgrey",
     lwd = 2, conf.lwd = 1, conf.lty = "dashed",
     shade = "lightgrey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object to plot; an object of class \code{"logitreg"}, usually
    the result of a call to \code{\link{logitreg}}.}
  \item{group}{The group to plot the logit model for. Can be one of the
    group labels or \code{"Combined"} to draw the individual logit
    models. Alternatively, and the default, is to use \code{"all"},
    which divides the plotting region into the required number of
    plotting regions and draws all the fitted curves.}
  \item{npred}{number of points at which the fitted curves are evaluated
    for plotting purposes.}
  \item{conf.int}{numeric; the confidence interval required.}
  \item{conf.type}{character; how should the confidence interval be
    drawn. Default is not to draw the confidence interval.}
  \item{xlab, ylab}{character; the x and y axis labels.}
  \item{rug}{logical; should rug plots be drawn?}
  \item{ticksize}{The size of the tick marks used in rug plots.}
  \item{col}{The colour in which to draw the line representing the
    fitted values.}
  \item{ref.col}{The colour of the reference lines drawn at 0 and 1.}
  \item{lwd}{The line width in which to draw the line representing the
    fitted values.}
  \item{conf.lwd, conf.lty}{Line width and line type for the confidence
    interval. Only used if \code{conf.type = "lines"} is specified.}
  \item{shade}{The colour for the fill and border of the confidence
    interval if \code{conf.type = "polygon"} is specified.}
  \item{\dots}{arguments passed on to \code{plot}.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A plot on the current device.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Gavin L. Simpson}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{logitreg}} for an example, \code{\link{roc}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
