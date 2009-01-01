\name{panel.Stratiplot}
\alias{panel.Stratiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Panel function for stratigraphic diagrams}
\description{
  A Lattice panel function for drawing individuals panels on
  stratigraphic diagrams using a range of plot types commonly used
  within palaeoecology.
}
\usage{
panel.Stratiplot(x, y,
                 type = "l",
                 col,
                 pch = plot.symbol$pch,
                 cex = plot.symbol$cex,
                 col.line = plot.line$col,
                 col.symbol = plot.symbol$col,
                 col.refline = ref.line$col,
                 col.smooth = "red",
                 col.poly = plot.line$col,
                 lty = plot.line$lty,
                 lwd = plot.line$lwd,
                 lty.smooth = plot.line$lty,
                 lwd.smooth = 2,
                 fill = plot.symbol$fill,
                 \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{variables defining the contents of the panel.}
  \item{type}{character vector consisting of one or more of the
    following: \code{"l"}, \code{"p"}, \code{"o"}, \code{"b"}, \code{"h"},
    \code{"g"}, \code{"smooth"}, and \code{"poly"}. It \code{type} has more
    than one element, the effects of each component are combined, though
    note that some elements will over-plot, and hence obscure, earlier
    elements.
  
    For \code{type}s \code{"l"}, \code{"p"}, \code{"o"}, \code{"b"} and
    \code{"g"} the standard Lattice interpretation is observed. See
    \code{\link{panel.xyplot}} for further details. Note that \code{type
      "b"} is the same as \code{type "o"}.

    \code{"g"} adds a reference grid using \code{panel.grid} in the
    background.
    
    For \code{"h"}, histogram-like bars are plotted from the
    \strong{y-axis}, not from the x-axis with \code{\link{plot}} and
    \code{\link{panel.loess}}.
  
    For \code{"smooth"} a loess fit is added to each panel using
    \code{\link{panel.Loess}}.
  
    For \code{"poly"}, a shaded polygon, or silhouette, is drawn for each
    panel.}
  \item{col, col.line, col.symbol, col.poly, col.refline,
    col.smooth}{colour parameters. For all but \code{col.smooth},
    default colours are obtained from \code{plot.symbol} and
    \code{plot.line} using
    \code{\link{trellis.par.get}}. \code{col.refline} controls the
    colour of the reference line drawn at value 0 on the x-axis.}
  \item{pch, cex, lty, lwd, fill}{other graphical parameters, defaults
    for which are obtained from \code{plot.symbol} and \code{plot.line}
    using \code{\link{trellis.par.get}}.}
  \item{lty.smooth}{line type for the loess smoother. The default is
    obtained from \code{plot.line} using \code{\link{trellis.par.get}}.}
  \item{lwd.smooth}{The line width for the loess smoother.}
  \item{\dots}{extra arguments passed on to the underlying panel
    functions; \code{\link{panel.points}}, \code{\link{panel.lines}},
    \code{\link{panel.segments}}, \code{\link{panel.polygon}},
    \code{\link{panel.Loess}} and \code{\link{panel.refline}}.}
}
\details{
  Creates stratigraphic scatter plots of \code{x} and \code{y}, with
  various modifications possible via the type argument.

  Note that all the arguments controlling the display can be supplied
  directly to a high-level call of the function \code{\link{Stratiplot}}.
}
\author{Gavin L. Simpson}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{Stratiplot}}, \code{\link{panel.Loess}},
    \code{\link{panel.xyplot}}.}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
