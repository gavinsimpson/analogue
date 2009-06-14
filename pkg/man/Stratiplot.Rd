\name{Stratiplot}
\alias{Stratiplot}
\alias{Stratiplot.default}
\alias{Stratiplot.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Palaeoecological stratigraphic diagrams}
\description{
  Draws palaeoecological stratigraphic diagrams of one or more variables
  as a function of depth/age, with the time dimension flowing from the
  bottom to the top of the y-axis, using the \pkg{Lattice} graphics
  package.
}
\usage{
Stratiplot(x, \dots)

\method{Stratiplot}{default}(x, y, type = "l", ylab = "", xlab = "",
           pages = 1, \dots)

\method{Stratiplot}{formula}(formula, data, subset, na.action, type = "l",
           ylab = "", xlab = "", pages = 1, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix-like object; the variables to be plotted.}
  \item{y}{numeric vector of depths/ages corresponding to rows in
    \code{x}. Length of \code{y} must be the same as \code{nrow(x)} or
    exactly equal to \code{nrow(x) / ncol(x)}. See Details.}
  \item{formula}{an object of class \code{"formula"} (or one that can be
    coerced to that class): a symbolic description of the model to be
    fitted. The details of plot specification are given under
    \sQuote{Details}.}
  \item{type}{character; The type of plotting. Can be a vector. Note
    that not all Lattice \sQuote{type}s are supported and some new types
    are allowed. See \code{\link{panel.Stratiplot}} for further
    details.}
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables to plot. If not found in \code{data}, the variables
    are taken from \code{environment(formula)}, typically the
    environment from which \code{Stratiplot} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be  used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s. The default is set by the \code{na.action}
    setting of \code{\link{options}}, and is \code{na.fail} if that is
    unset. The \sQuote{factory-fresh} default is \code{na.omit}. Another
    possible value is \code{NULL}, no action. Value \code{na.exclude} can be
    useful.}
  \item{ylab, xlab}{the x- and y-axis labels.}
  \item{pages}{numeric; the number of pages to draw the plot over. May
    be useful for data sets with many species.}
  \item{\dots}{additional arguments passed to
    \code{\link{panel.Stratiplot}} and the underlying
    \code{\link[lattice]{xyplot}} function.}
}
\details{
  Currently the function is designed for relative abundance data
  only. Greater flexibility is planned in a future version.

  Plots can be specified symbolically using a formula. A typical model
  has the form \code{Y ~ variables}, where \code{Y} is either the core
  depths or sample ages/dates (to be plotted on the y-axis) and
  \code{variables} is a series of terms which specifies the variables to
  plot against \code{Y}. Terms should be specified with the form
  \code{var1 + var2 + var3} to plot only those variables. Other,
  standard, notation for formulae apply, such as model formulae used in
  \code{\link{lm}}.

  Note that \code{formula} is \strong{not} passed on to
  \code{\link[lattice]{xyplot}}. Instead, the formula is parsed and evaluated
  within \code{Stratiplot} and an appropriate data structure formed to
  facilitate plotting via \code{\link[lattice]{xyplot}}. As such, the special
  features of \pkg{Lattice} formulae cannot be used.
}
\value{
  Returns a Lattice plot object of class \code{"trellis"}.
}
\author{Gavin L. Simpson.}
\note{
  The function currently doesn't know about ages/dates and will
  interpret these as \sQuote{depths} instead. This will be fixed in a
  future version.
}
\seealso{\code{\link[lattice]{xyplot}},
\code{\link{panel.Stratiplot}}, \code{\link{panel.Loess}}.}
\examples{
data(V12.122)
depths <- as.numeric(rownames(V12.122))

(plt <- Stratiplot(depths ~ O.univ + G.ruber + G.tenel + G.pacR,
                   data = V12.122,
                   type = c("h","l","g","smooth")))

(plt2 <- Stratiplot(depths ~ O.univ + G.ruber + G.tenel + G.pacR,
                    data = V12.122,
                    type = c("poly","g")))

(plt3 <- Stratiplot(depths ~ O.univ + G.ruber + G.tenel + G.pacR,
                    data = V12.122,
                    type = c("h")))

spp.want <- c("O.univ","G.ruber","G.tenel","G.pacR")
(plt4 <- Stratiplot(V12.122[, spp.want], y = depths,
                    type = c("poly", "g")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
