\name{bootstrapObject}
\alias{bootstrapObject}
\title{Bootstrap estimation and errors }
\description{
  Objects of class \code{\link{bootstrap.mat}} are a complex containing
  many sub-components. This object is described here in more detail.
}
\details{
  
  A large object is returned with some or all of the following depending
  on whether \code{newdata} and \code{newenv} are supplied or not.
  \describe{
    \item{\code{observed}: }{vector of observed environmental values.}
    \item{\code{model}: }{a list containing the apparent or non-bootstrapped
      estimates for the training set. With the following components:
      \describe{
	\item{\code{estimated}: }{estimated values for \code{"y"}, the
	  environment.}
	\item{\code{residuals}: }{model residuals.}
	\item{\code{r.squared}: }{Apparent \eqn{R^2} between observed and
	  estimated values  of \code{"y"}.}
	\item{\code{avg.bias}: }{Average bias of the model residuals.}
	\item{\code{max.bias}: }{Maximum bias of the model residuals.}
	\item{\code{rmse}: }{Apparent error (RMSE) for the model.}
	\item{\code{k}: }{numeric; indicating the size of model used in
	  estimates and predictions.}
      }
    }
    \item{\code{bootstrap}: }{a list containing the bootstrap estimates for the
      training set. With the following components:
      \describe{
	\item{\code{estimated}: }{Bootstrap estimates for \code{"y"}.}
	\item{\code{residuals}: }{Bootstrap residuals for \code{"y"}.}
	\item{\code{r.squared}: }{Bootstrap derived \eqn{R^2} between observed
	  and estimated values of \code{"y"}.}
	\item{\code{avg.bias}: }{Average bias of the bootstrap derived model
	  residuals.}
	\item{\code{max.bias}: }{Maximum bias of the bootstrap derived model
	  residuals.}
	\item{\code{rmsep}: }{Bootstrap derived RMSEP for the model.}
	\item{\code{s1}: }{Bootstrap derived S1 error component for the
	  model.}
	\item{\code{s2}: }{Bootstrap derived S2 error component for the
	  model.} 
	\item{\code{k}: }{numeric; indicating the size of model used in
	  estimates and predictions.}
      }
    }
    \item{\code{sample.errors}: }{a list containing the bootstrap-derived sample
      specific errors for the training set. With the following
      components:
      \describe{
	\item{\code{rmsep}: }{Bootstrap derived RMSEP for the training set
	  samples.} 
	\item{\code{s1}: }{Bootstrap derived S1 error component for training
	  set samples.}
	\item{\code{s2}: }{Bootstrap derived S2 error component for training
	  set samples.}
      }
    }
    \item{\code{weighted}: }{logical; whether the weighted mean was used instead of
      the mean of the environment for \emph{k}-closest analogues.}
    \item{\code{auto}: }{logical; whether \code{"k"} was choosen automatically or
      user-selected.}
    \item{\code{n.boot}: }{numeric; the number of bootstrap samples taken.}
    \item{\code{call}: }{the matched call.}
    \item{\code{type}: }{model type.}
    \item{\code{predictions}: }{a list containing the apparent and
      bootstrap-derived estimates for the new data, with the following
      components:
      \describe{
	\item{\code{observed}: }{the observed values for the new samples ---
	  only if \code{newenv} is provided.}
	\item{\code{model}: }{a list containing the apparent or
	  non-bootstrapped estimates for the new samples. A list with the
	  same components as \code{model}, above.}
	\item{\code{bootstrap}: }{a list containing the bootstrap estimates
	  for the new samples, with some or all of the same components as
	  \code{bootstrap}, above.}
	\item{\code{sample.errors}: }{a list containing the bootstrap-derived
	  sample specific errors for the new samples, with some or all of
	  the same components as \code{sample.errors}, above.}
      }
    }
  }
}
\author{Gavin L. Simpson }
\seealso{\code{mat}, \code{plot.mat}, \code{summary.bootstrap.mat},
  \code{\link{residuals}}}
\keyword{multivariate}
