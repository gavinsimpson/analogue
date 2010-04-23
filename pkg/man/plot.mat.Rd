\name{plot.mat}
\alias{plot.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot diagnostics for a mat object}
\description{
  Five plots (selectable by \code{which}) are currently available: a
  plot of estimated against observed values, a plot of residuals against
  estimated values, and screeplots of the apparent RMSE, average bias
  and maximum bias for MAT models of size \eqn{k}, where \eqn{k = 1,
    \dots, n}. By default, the first three and \sQuote{5} are provided.
}
\usage{
\method{plot}{mat}(x,
     which = c(1:3, 5),
     weighted = FALSE,
     k,
     caption = c("Inferred vs Observed", "Residuals vs Fitted",
                 "Leave-one-out errors", "Average bias",
                 "Maximum bias"),
     max.bias = TRUE,
     n.bias = 10,
     restrict = 20,
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
  \item{x}{an object of class \code{"mat"}.}
  \item{which}{which aspects of the \code{"mat"} object to plot if a
    subset of the plots is required, specify a subset of the numbers
    \code{1:5}.}
  \item{weighted}{logical; should the analysis use weighted mean of env
    data of analogues as fitted/estimated values?}
  \item{k}{numeric; the number of analogues to use. If missing \code{k}
    is chosen automatically as the \code{k} that achieves lowest RMSE.}
  \item{caption}{captions to appear above the plots.}
  \item{max.bias}{logical, should max bias lines be added to residuals.}
  \item{n.bias}{numeric, number of sections to calculate maximum bias
    for.}
  \item{restrict}{logical; restrict comparison of k-closest model to
    \eqn{k \le}{k <=} \code{restrict}.}
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
## Imbrie and Kipp example
## load the example data
data(ImbrieKipp)
data(SumSST)
data(V12.122)

## merge training and test set on columns
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- dat[[1]] / 100
V12.122 <- dat[[2]] / 100

## MAT
ik.mat <- mat(ImbrieKipp, SumSST, method = "chord")
ik.mat
summary(ik.mat)

## summary plot of MAT model
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(ik.mat)
layout(1)

}
\keyword{hplot}
\keyword{methods}
