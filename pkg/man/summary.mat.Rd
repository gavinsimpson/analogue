\name{summary.mat}
\alias{summary.mat}
\alias{print.summary.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarising Modern Analogue Technique models}
\description{
  \code{\link{summary}} method for class \code{"mat"}.
}
\usage{
\method{summary}{mat}(object, k = 10,
        digits = min(2, getOption("digits") - 4), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"cma"}, usually the
    result of a call to \code{\link{cma}}.}
  \item{k}{numeric; maximum modern analogues to use to summarise model
    fits.}
  \item{digits}{numeric; the number of significant digits with which to
    format results.}
  \item{\dots}{arguments passed to or from other methods.}
}
\value{
  A list with the components below. The number of analogues used,
  \emph{k} is returned as attribute \code{"k"}.

  \item{summ}{a data.frame containing the model fits for training set
    samples. See notes below.}
  \item{tbl}{matrix of summary statistics for an un-weighted model.}
  \item{tbl.W}{matrix of summary statistics for a weighted model.}
  \item{call}{the matched function call}
  \item{quantiles}{the quantiles of the distribution of pairwise
    dissimilarities for the training set, for \code{"probs = c(0.01,
      0.02, 0.05, 0.1, 0.2)"}.}
}
\author{Gavin L. Simpson}
\note{
  The returned component \code{"summ"} contains the following:

  \describe{
    \item{Obs:}{the observed responses for the training set samples.}
    \item{Est:}{the fitted values of the response for training set
      samples based on the average of \emph{k}-closest analogues.}
    \item{Resi:}{the residuals of the fitted model based on the average
      of \emph{k}-closest analogues.}
    \item{W.Est:}{the fitted values of the response for training set
      samples based on the weighted average of \emph{k}-closest
      analogues.}
    \item{W.Resi:}{the residuals of the fitted model based on the
      weighted average of \emph{k}-closest analogues.}
    \item{minDC:}{dissimilarity of closest analogue in training set for
      each training set sample.}
    \item{minResi:}{smallest residual for an un-weighted model of size
      \code{"k"}.}
    \item{k:}{size of model leading to minimal residual,
      \code{"minResi"}.}
    \item{minW.Resi:}{smallest residual for a weighted model of size
      \code{"k.W"}.}
    \item{k.W:}{size of model leading to minimal residual,
      \code{"minW.Resi"}.}
  }
}
\seealso{\code{\link{mat}}, \code{\link{summary}}.}
\examples{
\dontrun{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")
swap.mat

## model summary
summary(swap.mat)

## model summary - evaluating models using k = 1, ..., 20
## analogues instead of the default, 10.
summary(swap.mat, k = 20)
}
}
\keyword{methods}% at least one, from doc/KEYWORDS
