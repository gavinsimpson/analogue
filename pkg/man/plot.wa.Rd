\name{plot.wa}
\alias{plot.wa}
\title{Plot diagnostics for a weighted averaging model}
\description{
  Two plots (selectable by \code{which}) are currently available: a
  plot of estimated against observed values, a plot of residuals against
  estimated values.
}
\usage{
\method{plot}{wa}(x,
     which = 1:2,
     caption = c("Inferred vs Observed", "Residuals vs Fitted"),
     max.bias = TRUE,
     n.bias = 10,
     sub.caption = NULL,
     main = "",
     ask = prod(par("mfcol")) < length(which) &&
                                  dev.interactive(),
     \dots,
     panel = if (add.smooth) panel.smooth else points,
     add.smooth = getOption("add.smooth"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"wa"}.}
  \item{which}{which aspects of the \code{"wa"} object to plot if a
    subset of the plots is required, specify a subset of the numbers
    \code{1:2}.}
  \item{caption}{captions to appear above the plots.}
  \item{max.bias}{logical, should max bias lines be added to residuals.}
  \item{n.bias}{numeric, number of sections to calculate maximum bias
    for.}
  \item{sub.caption}{common title-above figures if there are multiple;
    used as \sQuote{sub} (s.\sQuote{title}) otherwise.  If \code{NULL},
    as by default, a possibly shortened version of
    \code{deparse(x$call)} is used.}
  \item{main}{title to each plot-in addition to the above
    \code{caption}.}
  \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, see \code{par(ask=.)}.}
  \item{\dots}{graphical arguments passed to other graphics functions.}
  \item{panel}{panel function.  The useful alternative to
    \code{points}, \code{panel.smooth}, can be chosen by
    \code{add.smooth = TRUE}.}
  \item{add.smooth}{logical indicating if a smoother should be added to
    fitted \& residuals plots; see also \code{panel} above.}
}
\details{
  This plotting function is modelled closely on \code{\link{plot.lm}}
  and many of the conventions and defaults for that function are
  replicated here.
  
  \code{sub.caption} - by default the function call - is shown as a
  subtitle (under the x-axis title) on each plot when plots are on
  separate pages, or as a subtitle in the outer margin (if any) when
  there are multiple plots per page.
}
\value{
  One or more plots, drawn on the current device.
}
%\references{FIXME: one or two references where these types of plot are
%  discussed in the literature.}
\author{Gavin L. Simpson. Code borrows heavily from \code{\link{plot.lm}}.}
\seealso{\code{\link{mat}}}
\examples{
## see full example in ?wa

}
\keyword{hplot}
\keyword{methods}
