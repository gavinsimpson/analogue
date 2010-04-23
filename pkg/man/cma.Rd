\name{cma}
\alias{cma}
\alias{cma.default}
\alias{cma.analog}
\alias{print.cma}
\alias{plot.cma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Close modern analogues}
\description{
  Extracts and formats close modern analogue samples from a modern
  reference set that are closer than a defined cut off threshold.
}
\usage{
cma(object, \dots)

\method{cma}{default}(object, \dots)

\method{cma}{analog}(object, cutoff, prob = c(0.01, 0.025, 0.05), \dots)

\method{plot}{cma}(x, method = c("overplot", "jitter", "stack"),
   jitter = 0.1, vertical = FALSE,
   draw.quant = TRUE, xlab = NULL, ylab = "",
   main = "", cex.axis = NULL, \dots,
   col.quant = "red", lty.quant= "dashed")
}
\arguments{
  \item{object}{an object for which close modern analogues are to be
    returned. Currently only for objects of class \code{\link{analog}}.}
  \item{cutoff}{numeric; critical value determining level below which
    samples from the modern reference set are defined as close modern
    analogues. May be missing, in which case the 2.5\% quantile of the
    training set dissimilarities is used unless \code{object$train} is
    \code{NULL}, in which case \code{"cutoff"} must be supplied.}
  \item{prob}{numeric vector of probabilities with values in [0,1], for
    which quantiles of the distribution of training set dissimilarities
    will be calculated. See \code{\link{quantile}}.}
  \item{\dots}{arguments to be passed to other \code{\link{cma}}
    methods or additional arguments passed to \code{\link{stripchart}}.}
  \item{x}{an object of class \code{"cma"}.}
  \item{method}{the method to be used to separate coincident points. The
    default method\cr\code{"overplot"} causes such points to be
    overplotted, but it is also possible to specify \code{"jitter"} to 
    jitter the points, or \code{"stack"} have coincident points
    stacked.  The last method only makes sense for very granular data.}
  \item{jitter}{when \code{method="jitter"} is used, \code{jitter} gives
    the amount of jittering applied.}
  \item{vertical}{when vertical is \code{TRUE} the plots are drawn
    vertically rather than the default horizontal.}
  \item{draw.quant}{logical; should the quantiles be drawn on the stripchart?}
  \item{xlab,ylab,main}{Graphical parameters}
  \item{cex.axis}{The magnification to be used for axis annotation
    relative to the current setting of \code{cex}. See
    \code{\link{par}}.}
  \item{col.quant,lty.quant}{colour and line type in which to drawn the
    quantile lines.}
}
\details{
  The plot method is simply a wrapper to \code{\link{stripchart}}.
}
\value{
  For the plot method, a plot on the current device. Invisibly the
  plotted data are returned; see Note for further details.

  A list of class \code{"cma"} with the following components:
  \item{close}{a named list of named vectors of close modern analogues
    and their dissimilarities. The names of the list components are the
    names of the fossil samples. The named vector in each
    component of \code{close} is the distances for the close modern
    analogues from the training set that are as close as \code{cutoff},
    or closer, to the fossil sample.}
  \item{call}{the matched call.}
  \item{cutoff}{the cutoff threshold used to define close modern
    analogues.}
  \item{quant}{numeric vector of the requested quantiles.}
  \item{probs}{the probabilities of the requested quantiles.}
  \item{method}{character; the dissimilarity coefficient used}
  \item{n.analogs}{numeric vector of the number of analogues per fossil
    sample.}
}
\note{
  Currently, only objects of class \code{\link{analog}} are supported.
  
  The plot method invisibly returns a list with the following
  components:
  \itemize{
    \item{\code{distances}}{a vector of stacked distances extracted from
      \code{object}.}
    \item{\code{groups}}{a factor listing the fossil sample for which the
      distances are the distances to the close modern analogues for the
      training set.}
  }
}
\references{
  
  Flower, R.J., Juggins, S. and Battarbee, R.W. (1997) Matching diatom
  assemblages in lake sediment cores and modern surface sediment
  samples: the implications for lake conservation and restoration with
  special reference to acidified systems. \emph{Hydrobiologia}
  \strong{344}; 27--40.

  Simpson, G.L., Shilland, E.M., Winterbottom, J. M. and Keay, J. (2005)
  Defining reference conditions for acidified waters using a modern
  analogue approach. \emph{Environmental Pollution} \strong{137};
  119--133.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{analog}}, \code{\link{stripchart}}, or
  \code{\link{boxplot}} for an alternative representation.}
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

## analog matching between SWAP and RLGH reference samples
(ik.ana <- analog(ImbrieKipp, V12.122, method = "chord"))

## close modern analogues --- no analogues are "close"
(ik.cma <- cma(ik.ana, cutoff = 0.91))
summary(ik.cma)

## plot the results
plot(ik.cma)

}
\keyword{methods}
\keyword{manip}
\keyword{hplot}