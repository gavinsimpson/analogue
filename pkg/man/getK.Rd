\name{getK}
\alias{getK}
\alias{getK.default}
\alias{getK.mat}
\alias{getK.bootstrap.mat}
\alias{getK.predict.mat}
\alias{setK<-}
\alias{setK<-.default}
\alias{setK<-.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract and set the number of analogues}
\description{
  An extractor function to access the number of analogues used in
  particular models. The stored value of \eqn{k} can be updated using
  \code{setK}.
}
\usage{
getK(object, ...)

\method{getK}{mat}(object, weighted = FALSE, ...)

\method{getK}{bootstrap.mat}(object, which = c("bootstrap", "model"),
     prediction = FALSE, ...)

\method{getK}{predict.mat}(object, which = c("model", "bootstrap"),
     ...)

setK(object, weighted = FALSE) <- value

\method{setK}{mat}(object, weighted = FALSE) <- value

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object; currently only for objects of class
    \code{\link{mat}} and class \code{\link{bootstrap.mat}}.}
  \item{weighted}{logical; extract/set number of analogues for a
    weighted or un-weighted model?}
  \item{which}{character; which \emph{k} should be extracted, the one
    from the model or the one from the bootstrap results?}
  \item{prediction}{logical; should the extracted \emph{k} be the one
    that is minimum for the test set (\code{newdata}) or the model
    (\code{object}).}
  \item{\dots}{further arguments to other methods.}
  \item{value}{integer; replacement value for \eqn{k}.}
}
\details{
  \code{getK} is a generic accessor function, and \code{setK<-} is a generic
  replacement function.

  Objects of class \code{\link{bootstrap.mat}} contain several different
  \code{k}'s. If no predictions are performed, there will be two
  \code{k}'s, one for the model and one from bootstrapping the
  model. Where predictions are performed \strong{with} \code{newenv}
  supplied, in addition to the \code{k}'s above, there will be two
  \code{k}' for the predictions, one for the model-based and one for the
  bootstrap-based predictions. To select \code{k} for the predictions,
  use \code{prediction = TRUE}. Argument \code{which} determines whether
  the model-based or the bootstrap-based \code{k} is returned.
}
\value{
  For \code{getK}, an integer value that is the number of analogues stored
  for use. The returned object has attributes \dQuote{auto} and
  \dQuote{weighted}. \dQuote{auto} refers to whether the extracted value
  of \eqn{k} was set automatically (\code{TRUE}) or by the user
  (\code{FALSE}). \dQuote{weighted} states if the returned value is for
  a \code{weighted} analysis or an un-\code{weighted} analysis (\code{FALSE}).

  For \code{setK<-}, the updated object.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{mat}}}
\examples{
## continue the example from join
example(join)

## fit a MAT model
ik.mat <- mat(ImbrieKipp, SumSST, method = "SQchord")

## How many analogues gives lowest RMSE?
getK(ik.mat)
## note that this value was chosen automatically

## Now set k to be 10
setK(ik.mat) <- 10

## check
getK(ik.mat)

}
\keyword{utilities}
\keyword{manip}
