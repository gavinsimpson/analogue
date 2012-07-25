\name{panel.Stratiplot}
\alias{panel.Stratiplot}

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
                 lwd.h = 3,
                 fill = plot.symbol$fill,
                 zones = NULL,
                 col.zones = plot.line$col,
                 lty.zones = plot.line$lty,
                 lwd.zones = plot.line$lwd,
                 gridh = -1, gridv = -1,
                 \dots)
}

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
    \code{\link[lattice]{panel.xyplot}} for further details. Note that
    \code{type "b"} is the same as \code{type "o"}.

    \code{"g"} adds a reference grid using \code{panel.grid} in the
    background.
    
    For \code{"h"}, histogram-like bars are plotted from the
    \strong{y-axis}, not from the x-axis with \code{\link{plot}} and
    \code{\link[lattice]{panel.loess}}.
  
    For \code{"smooth"} a loess fit is added to each panel using
    \code{\link{panel.Loess}}.
  
    For \code{"poly"}, a shaded polygon, or silhouette, is drawn for each
    panel.}
  \item{col, col.line, col.symbol, col.poly, col.refline,
    col.smooth, col.zones}{colour parameters. For all but \code{col.smooth},
    default colours are obtained from \code{plot.symbol} and
    \code{plot.line} using
    \code{\link[lattice]{trellis.par.get}}. \code{col.refline} controls the
    colour of the reference line drawn at value 0 on the x-axis, as well
    as the colour of the grid lines if drawn.}
  \item{pch, cex, lty, lwd, fill}{other graphical parameters, defaults
    for which are obtained from \code{plot.symbol} and \code{plot.line}
    using \code{\link[lattice]{trellis.par.get}}.}
  \item{lty.smooth}{line type for the loess smoother. The default is
    obtained from \code{plot.line} using
    \code{\link[lattice]{trellis.par.get}}.}
  \item{lwd.smooth,lwd.h}{The line width for the loess smoother and
    histogram-like bars respectively.}
  \item{zones}{numeric; vector of zone boundary positions on scale of
    the depth/time (y-)axis.}
  \item{lty.zones, lwd.zones}{line type and width for the zone
    markers. The defaults are obtained from \code{plot.line}.}
  \item{gridh, gridv}{numeric arguments corresponding to \code{h} and
    \code{v} of \code{\link{panel.grid}}, which control the number of
    grid lines drawn.}
  \item{\dots}{extra arguments passed on to the underlying panel
    functions; \code{\link[lattice]{panel.points}},
    \code{\link[lattice]{panel.lines}},
    \code{\link[lattice]{panel.segments}},
    \code{\link[lattice]{panel.polygon}}, 
    \code{\link{panel.Loess}} and \code{\link[lattice]{panel.refline}}.}
}
\details{
  Creates stratigraphic scatter plots of \code{x} and \code{y}, with
  various modifications possible via the type argument.

  Zones can be drawn on the panels by supplying the numeric vector of
  zone boundaries as argument \code{zones}. The panel function will then
  draw horizontal lines across the panels at the desired y-axis
  locations. Note that the panel function does \strong{not} attempt to
  identify the zone boundaries automatically; these must be determined
  via a chronological (constrained) cluster analysis function or
  similar.

  Note that all the arguments controlling the display can be supplied
  directly to a high-level call of the function \code{\link{Stratiplot}}.
}

\note{
  Histogram-like bars (\code{type = "h"}) are drawn with \code{lineend =
  "butt"} to improve their appearance. This can not be changed by the
  user and you can't include that grid parameter in any call to
  \code{panel.Stratiplot} that uses \code{type = "h"}.
}

\author{Gavin L. Simpson}

\seealso{\code{\link{Stratiplot}}, \code{\link{panel.Loess}},
    \code{\link[lattice]{panel.xyplot}}.}
\keyword{hplot}
