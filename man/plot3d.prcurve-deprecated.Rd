\name{plot3d.prcurve-deprecated}
\alias{plot3d.prcurve}

\title{Interactive 3D plot of a principal curve in principal coordinate
  space}
\description{
  Draws a 3D plot of the principal curve in principal coordinate space
  using the \pkg{rgl} package and functions from \pkg{vegan}.
}
\usage{
plot3d.prcurve(x, choices = 1:3, display = "sites", scaling = 0,
               lcol = "darkorange", lwd = 2, decorate = TRUE, xlab = NULL,
               ylab = NULL, zlab = NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object of class \code{"prcurve"} resulting from a call to
    \code{\link{prcurve}}.
  }
  \item{choices}{
    numeric vector of length 3; the ordination axes to plot.
  }
  \item{display}{
    character; which scores to display. See \code{\link{scores.rda}}.
  }
  \item{scaling}{
    numeric; the scaling to use for the scores. Default is no scaling.
  }
  \item{lcol, lwd}{
    The colour and width, respectively, for the principal curve.
  }
  \item{decorate}{
    logical; should the plot be decorated with bounding box, axes and
    labels?
  }
  \item{xlab, ylab, zlab, main}{
    Labels for the plot, to be drawn using \code{decorate3d}.
  }
  \item{\dots}{
    Arguments passed to other functions. In particular, argments are
    passed to \code{\link{ordirgl}}, \code{lines3d}, and
    \code{decorate3d}.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A plot is drawn on the active RGL device. If there is no active RGL
  device, one is opened upon plotting.
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
  \code{\link{ordirgl}}, \code{\link{prcurve}}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dynamic }
\keyword{ graphics }% __ONLY ONE__ keyword per line
