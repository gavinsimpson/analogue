\name{predict.mat}
\alias{predict.mat}
\alias{print.predict.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict method for Modern Analogue Technique models}
\description{
  Predicted values based on a MAT model object.
}
\usage{
\method{predict}{mat}(object, newdata, k, weighted = FALSE,
        bootstrap = FALSE, n.boot = 1000,
        probs = c(0.01, 0.025, 0.05, 0.1), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of \code{\link{mat}}.}
  \item{newdata}{data frame; required only if predictions for some new
    data are required. Mst have the same number of columns, in same
    order, as \code{x} in \code{\link{mat}}. See example below or
    \code{\link{join}} for how to do this. If \code{newdata} not
    provided, the fitted values are returned.}
  \item{k}{number of analogues to use. If missing, \code{k} is chosen
    automatically as the \code{k} that achieves lowest RMSE.}
  \item{weighted}{logical; should the analysis use the weighted mean of
    environmental data of analogues as predicted values?}
  \item{bootstrap}{logical; should bootstrap derived estimates and
    sample specific errors be calculated-ignored if \code{newdata} is
    missing.}
  \item{n.boot}{numeric; the number of bootstrap samples to take.}
  \item{probs}{numeric; vector of probabilities with values in [0,1].}
  \item{\dots}{arguments passed to of from other methods.}
}
\details{
  This function returns one or more of three sets of results depending
  on the supplied arguments:
  \describe{
    \item{Fitted values:}{the fitted values of the \code{\link{mat}}
      model are returned if \code{newdata} is missing.}
    \item{Predicted values:}{the predicted values for some new samples
      are returned if \code{newdata} is supplied. Summary model
      statistics and estimated values for the training set are also
      returned.}
    \item{Bootstrap predictions and standard errors:}{if \code{newdata}
      is supplied and \code{bootstrap = TRUE}, the predicted values for
      \code{newdata} plus bootstrap estimates and standard errors for the
      new samples and the training set are returned.}
  }
  The latter is simply a wrapper for \code{bootstrap(model, newdata,
    \dots)} - see \code{\link{bootstrap.mat}}.
}
\value{
  A object of class \code{predict.mat} is returned if \code{newdata} is
  supplied, otherwise an object of \code{\link{fitted.mat}} is
  returned. If \code{bootstrap = FALSE} then not all components will be
  returned.

  \item{observed}{vector of observed environmental values.}
  \item{model}{a list containing the model or non-bootstrapped
    estimates for the training set. With the following components:
    \itemize{
      \item{\code{estimated}}{estimated values for \code{"y"}, the
	environment.}
      \item{\code{residuals}}{model residuals.}
      \item{\code{r.squared}}{Model \eqn{R^2} between observed and
	estimated values  of \code{"y"}.}
      \item{\code{avg.bias}}{Average bias of the model residuals.}
      \item{\code{max.bias}}{Maximum bias of the model residuals.}
      \item{\code{rmsep}}{Model error (RMSEP).}
      \item{\code{k}}{numeric; indicating the size of model used in
	estimates and predictions.}
    }
  }
  \item{bootstrap}{a list containing the bootstrap estimates for the
    training set. With the following components:
    \itemize{
      \item{\code{estimated}}{Bootstrap estimates for \code{"y"}.}
      \item{\code{residuals}}{Bootstrap residuals for \code{"y"}.}
      \item{\code{r.squared}}{Bootstrap derived \eqn{R^2} between observed
	and estimated values of \code{"y"}.}
      \item{\code{avg.bias}}{Average bias of the bootstrap derived model
	residuals.}
      \item{\code{max.bias}}{Maximum bias of the bootstrap derived model
	residuals.}
      \item{\code{rmsep}}{Bootstrap derived RMSEP for the model.}
      \item{\code{s1}}{Bootstrap derived S1 error component for the
	model.}
      \item{\code{s2}}{Bootstrap derived S2 error component for the
	model.} 
      \item{\code{k}}{numeric; indicating the size of model used in
	estimates and predictions.}
    }
  }
  \item{sample.errors}{a list containing the bootstrap-derived sample
    specific errors for the training set. With the following components:
    \itemize{
      \item{\code{rmsep}}{Bootstrap derived RMSEP for the training set
	samples.} 
      \item{\code{s1}}{Bootstrap derived S1 error component for training
	set samples.}
      \item{\code{s2}}{Bootstrap derived S2 error component for training
	set samples.}
    }
  }
  \item{weighted}{logical; whether the weighted mean was used instead of
    the mean of the environment for \emph{k}-closest analogues.}
  \item{auto}{logical; whether \code{"k"} was choosen automatically or
    user-selected.}
  \item{n.boot}{numeric; the number of bootstrap samples taken.}
  \item{predictions}{a list containing the model and
    bootstrap-derived estimates for the new data, with the following
    components:
    \itemize{
      \item{\code{observed}}{the observed values for the new samples ---
	only if \code{newenv} is provided.}
      \item{\code{model}}{a list containing the model or
	non-bootstrapped estimates for the new samples. A list with the
	same components as \code{model}, above.
      }
      \item{\code{bootstrap}}{a list containing the bootstrap estimates
	for the new samples, with some or all of the same components as
	\code{bootstrap}, above.}
      \item{\code{sample.errors}}{a list containing the bootstrap-derived
	sample specific errors for the new samples, with some or all of
	the same components as \code{sample.errors}, above.}
    }
  }
}

\references{
  
  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C. and ter Braak,
  C.J.F. (1990). Diatoms and pH reconstruction. \emph{Philosophical
    Transactions of the Royal Society of London; Series B}, \strong{327};
  263--278.
}

\author{Gavin L. Simpson}
\seealso{\code{\link{mat}}, \code{\link{bootstrap.mat}}}
\examples{
## continue the RLGH and SWAP example from ?join
example(join)

## fit the MAT model using the squared chord distance measure
swap.mat <- mat(swapdiat, swappH, method = "SQchord")

## predict for RLGH data
predict(swap.mat, rlgh)

}
\keyword{methods}% at least one, from doc/KEYWORDS
