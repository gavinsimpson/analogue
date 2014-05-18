\name{summary.bootstrap.mat}
\alias{summary.bootstrap.mat}
\alias{print.summary.bootstrap.mat}
%\name{summary.bootstrap}
%\alias{summary.bootstrap}
%\alias{print.summary.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarise bootstrap resampling for MAT models}
\description{
  \code{\link{summary}} method for class \code{"bootstrap.mat"}.
}
\usage{
\method{summary}{bootstrap.mat}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"bootstrap.mat"}, usually the
    result of a call to \code{\link{bootstrap.mat}}.}
  \item{\dots}{arguments passed to or from other methods.}
}
\value{
  A data frame with the following components:
  
  \item{observed}{vector of observed environmental values.}
  \item{model}{a list containing the apparent or non-bootstrapped
    estimates for the training set. With the following components:
    \describe{
      \item{\code{estimated}: }{estimated values for the response}
      \item{\code{residuals}: }{model residuals}
      \item{\code{r.squared}: }{Apparent \eqn{R^2} between observed and
	estimated values  of \code{y}}
      \item{\code{avg.bias}: }{Average bias of the model residuals}
      \item{\code{max.bias}: }{Maximum bias of the model residuals}
      \item{\code{rmse}: }{Apparent error (RMSE) for the model}
      \item{\code{k}: }{numeric; indicating the size of model used in
	estimates and predictions}
    }
  }
  \item{bootstrap}{a list containing the bootstrap estimates for the
    training set. With the following components:
    \describe{
      \item{\code{estimated}: }{Bootstrap estimates for the response}
      \item{\code{residuals}: }{Bootstrap residuals for the response}
      \item{\code{r.squared}: }{Bootstrap derived \eqn{R^2} between observed
	and estimated values of the response}
      \item{\code{avg.bias}: }{Average bias of the bootstrap derived model
	residuals}
      \item{\code{max.bias}: }{Maximum bias of the bootstrap derived model
	residuals}
      \item{\code{rmsep}: }{Bootstrap derived RMSEP for the model}
      \item{\code{s1}: }{Bootstrap derived S1 error component for the
	model}
      \item{\code{s2}: }{Bootstrap derived S2 error component for the
	model} 
      \item{\code{k}: }{numeric; indicating the size of model used in
	estimates and predictions}
    }
  }
  \item{sample.errors}{a list containing the bootstrap-derived sample
    specific errors for the training set. With the following components:
    \describe{
      \item{\code{rmsep}: }{Bootstrap derived RMSEP for the training set
	samples} 
      \item{\code{s1}: }{Bootstrap derived S1 error component for training
	set samples}
      \item{\code{s2}: }{Bootstrap derived S2 error component for training
	set samples}
    }
  }
  \item{weighted}{logical; whether the weighted mean was used instead of
    the mean of the environment for \emph{k}-closest analogues}
  \item{auto}{logical; whether \code{k} was choosen automatically or
    user-selected}
  \item{n.boot}{numeric; the number of bootstrap samples taken}
  \item{call}{the matched call}
  \item{call}{model type}
  \item{predictions}{a list containing the apparent and
    bootstrap-derived estimates for the new data, with the following
    components:
    \describe{
      \item{\code{observed}: }{the observed values for the new samples ---
	only if \code{newenv} is provided}
      \item{\code{model}: }{a list containing the apparent or
	non-bootstrapped estimates for the new samples. A list with the
	same components as \code{apparent}, above}
      \item{\code{bootstrap}: }{a list containing the bootstrap estimates
	for the new samples, with some or all of the same components as
	\code{bootstrap}, above}
      \item{\code{sample.errors}: }{a list containing the bootstrap-derived
	sample specific errors for the new samples, with some or all of
	the same components as \code{sample.errors}, above}
    }
  }
}
\author{Gavin L. Simpson}
\seealso{\code{\link{bootstrap.mat}}, \code{\link{mat}},
  \code{\link{summary}}.}
\examples{
\dontrun{
## continue the RLGH example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## bootstrap training set
swap.boot <- bootstrap(swap.mat, k = 10, n.boot = 100)
swap.boot
summary(swap.boot)
}
}
\keyword{methods}% at least one, from doc/KEYWORDS
