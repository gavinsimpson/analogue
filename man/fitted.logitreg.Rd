\name{fitted.logitreg}
\alias{fitted.logitreg}

\title{Fitted values for the training set from logistic regression
  models}

\description{
  Extracts fitted values for training set samples from logistic
  regression models fitted to each group of samples that describe the
  probability two samples are analogues (from the same group) as a
  function of dissimilarity between the paired samples. 
}
\usage{
\method{fitted}{logitreg}(object, combined = FALSE, ...)
}

\arguments{
  \item{object}{an object of class \code{"logitreg"} resulting from a
    call to \code{\link{logitreg}}.}
  \item{combined}{logical; should the fitted values for the overall
    combined analysis be returned.}
  \item{\dots}{arguments passed to other methods.}
}

\value{
  If \code{combined == FALSE} (the default) then a matrix of fitted
  probabilities, where the rows are the training set samples and the
  columns the groupings, is returned. If \code{combined == TRUE}, then a
  list with components \code{"group"} and
  \code{"combined"}. \code{"group"} is a matrix of fitted probabilities
  as above. \code{"combined"} is a vector of fitted values for the
  entire set of pairwise comparisons considered.
}

\author{Gavin L. Simpson}

\seealso{
  See \code{\link{logitreg}} for example usage.
}

\keyword{methods}
