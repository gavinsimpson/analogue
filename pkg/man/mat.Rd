\name{mat}
\alias{mat}
\alias{mat.default}
\alias{mat.formula}
\alias{fitted.mat}
\alias{residuals.mat}
\alias{resid.mat}
\alias{print.residuals.mat}
\alias{print.mat}
\alias{print.fitted.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Modern Analogue Technique transfer function models}
\description{
  Modern Analogue Technique (MAT) transfer function models for
  palaeoecology. The fitted values are the, possibly weighted, averages
  of the environment for the \emph{k}-closest modern analogues. MAT is a
  \emph{k}-NN method.
}
\usage{
mat(x, \dots)

\method{mat}{default}(x, y,
    method = c("euclidean", "SQeuclidean", "chord", "SQchord",
               "bray", "chi.square", "SQchi.square",
               "information", "chi.distance", "manhattan",
               "kendall", "gower", "alt.gower", "mixed"),
    kmax, \dots)

\method{mat}{formula}(formula, data, subset, na.action,
    method = c("euclidean", "SQeuclidean", "chord", "SQchord",
               "bray", "chi.square", "SQchi.square",
               "information", "chi.distance", "manhattan",
               "kendall", "gower", "alt.gower", "mixed"),
    model = FALSE, \dots)

\method{fitted}{mat}(object, k, weighted = FALSE, \dots)

\method{residuals}{mat}(object, k, weighted = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data frame containing the training set data, usually
    species data.}
  \item{y}{a vector containing the response variable, usually
    environmental data to be predicted from \code{x}.}
  \item{formula}{a symbolic description of the model to be fit. The
    details of model specification are given below.}
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)}, typically the
    environment from which \code{wa} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  The default is set by the
    \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The "factory-fresh" default
    is \code{\link{na.omit}}.  Another possible value is \code{NULL}, no
    action. Value \code{\link{na.exclude}} can be useful.}
  \item{method}{a character string indicating the dissimilarity
    (distance) coefficient to be used to define modern analogues. See
    Details, below.}
  \item{model}{logical; If \code{TRUE} the model frame of the fit is
    returned.}
  \item{kmax}{numeric; limit the maximum number of analogues considered
    during fitting. By default, \code{kmax} is equal to \eqn{n - 1},
    where \eqn{n} is the number of sites. For large data sets this is
    just wasteful as we wouldn't expect to be averaging over the entire
    training set. \code{kmax} can be used to restrict the upper limit on
    the number of analogues considered.}
  \item{object}{an object of class \code{\link{mat}}.}
  \item{k}{numeric; the \emph{k}-closest analogue models' for which
    fitted values and residuals are returned. Overides the default stored in
    the object.}
  \item{weighted}{logical; should weighted averages be used instead of
    simple averages?}
  \item{\dots}{arguments can be passed to \code{\link{distance}} to
    provide additional optios required for some dissimilarities.}
}
\details{
  The Modern Analogue Technique (MAT) is perhaps the simplest of the
  transfer function models used in palaeoecology. An estimate of the
  environment, \eqn{x}, for the response for a fossil sample, \eqn{y},
  is the, possibly weighted, mean of that variable across the
  \emph{k}-closest modern analogues selected from a modern  training set
  of samples. If used, weights are the reciprocal of the dissimilarity
  between the fossil sample and each modern analogue.

  A typical model has the form \code{response ~ terms} where
  \code{response} is the (numeric) response data frame and \code{terms}
  is a series of terms which specifies a linear predictor for
  \code{response}. A typical form for \code{terms} is \code{.},
  which is shorthand for "all variables" in \code{data}. If \code{.} is
  used, \code{data} must also be provided. If specific species
  (variables) are required then \code{terms} should take the form
  \code{spp1 + spp2 + spp3}.

  Pairwise sample dissimilarity is defined by dissimilarity or
  distance coefficients. A variety of coefficients are supported --- see
  \code{\link{distance}} for details of the supported coefficients.

  \emph{k} is chosen by the user. The simplest choice for \emph{k} is to
  evaluate the RMSE of the difference between the predicted and observed
  values of the environmental variable of interest for the training set
  samples for a sequence of models with increasing \emph{k}. The number
  of analogues chosen is the value of \emph{k} that has lowest
  RMSE. However, it should be noted that this value is biased as the
  data used to build the model are also used to test the predictive
  power.

  An alternative approach is to employ an optimisation data set on which
  to evaluate the size of \eqn{k} that provides the lowest RMSEP. This
  may be impractical with smaller sample sizes.

  A third option is to bootstrap re-sample the training set many times. At
  each bootstrap sample, predictions for samples in the bootstrap test
  set can be made for \eqn{k = 1, ..., n}, where \eqn{n} is the
  number of samples in the training set. \eqn{k} can be chosen from the
  model with the lowest RMSEP. See function \code{\link{bootstrap.mat}} for
  further details on choosing \eqn{k}.

  The output from \code{\link{summary.mat}} can be used to choose
  \eqn{k} in the first case above. For predictions on an optimsation or
  test set see \code{\link{predict.mat}}. For bootstrap resampling of
  \code{\link{mat}} models, see \code{\link{bootstrap.mat}.}
  
  The fitted values are for the training set and are taken as the,
  possibly weighted, mean of the environmental variable in question
  across the \emph{k}-closest analogues. The fitted value for each
  sample does \bold{not} include a contribution from itself --- it is
  the closest analogue, having zero dissimilarity. This spurious
  distance is ignored and analogues are ordered in terms of the non-zero
  distances to other samples in the training set, with the
  \emph{k}-closest contributing to the fitted value.

  Typical usages for \code{residuals.mat} are:
  \preformatted{
    resid(object, k, weighted = FALSE, \dots)
  }
}
\value{
  \code{mat} returns an object of class \code{mat} with the following
  components: 
  \item{standard}{list; the model statistics based on simple
    averages of \emph{k}-closest analogues. See below.}
  \item{weighted}{list; the model statistics based on weighted  of
    \emph{k}-closest analogues. See below.}
  \item{Dij}{matrix of pairwise sample dissimilarities for the training
    set \code{x}.}
  \item{orig.x}{the original training set data.}
  \item{orig.y}{the original environmental data or response, \code{y}.}
  \item{call}{the matched function call.}
  \item{method}{the dissimilarity coefficient used.}

  If \code{model = TRUE} then additional components \code{"terms"} and
  \code{"model"} are returned containing the \code{\link{terms}} object
  and model frame used.
  
  \code{fitted.mat} returns a list with the following components:
  \itemize{
    \item{estimated}{numeric; a vector of fitted values.}
    \item{k}{numeric; this is the \emph{k}-closest analogue model with
      lowest apparent RMSE.}
    \item{weighted}{logical; are the fitted values the weighted averages
      of the environment for the \emph{k}-closest analogues. If
      \code{FALSE}, the fitted values are the average of the environment
      for the \emph{k}-closest analogues.}
  }
}
\note{
  The object returned by \code{mat} contains lists \code{"standard"} and
  \code{"weighted"} both with the following elements:
  \itemize{
    \item{\code{est}}{a matrix of estimated values for the training set
      samples for models using \eqn{k} analogues, where \eqn{k = 1,
	..., n}. \eqn{n} is the number of smaples in the training
      set. Rows contain the values for each model of size \eqn{k}, with
      colums containing the estimates for each training set sample.}
    \item{\code{resid}}{matrix; as for \code{"est"}, but containing the
      model residuals.}
    \item{\code{rmsep}}{vector; containing the leave-one-out root mean square
      error or prediction.}
    \item{\code{avg.bias}}{vector; contains the average bias (mean of
      residuals) for models using \emph{k} analogues, where \eqn{k = 1,
	..., n}. \eqn{n} is the number of smaples in the training set.}
    \item{\code{max.bias}}{vector; as for \code{"avg.bias"}, but
      containing the maximum bias statistics.}
    \item{\code{r.squared}}{vector; as for \code{"avg.bias"}, but
      containing the \eqn{R^2} statistics.}
    }
}
\references{
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.
  
  Overpeck, J.T., Webb III, T. and Prentice I.C. (1985) Quantitative
  interpretation of fossil pollen spectra: dissimilarity coefficients and
  the method of modern analogues. \emph{Quaternary Research} \strong{23},
  87--108.

  Prell, W.L. (1985) The stability of low-latitude sea-surface
  temperatures: an evaluation of the CLIMAP reconstruction with emphasis
  on the positive SST anomalies, Report TR 025. U.S. Department of
  Energy, Washington, D.C.
  
  Sawada, M., Viau, A.E., Vettoretti, G., Peltier, W.R. and Gajewski,
  K. (2004) Comparison of North-American pollen-based temperature and
  global lake-status with CCCma AGCM2 output at 6 ka. \emph{Quaternary
    Science Reviews} \strong{23}, 87--108.

}
\author{Gavin L. Simpson}
\seealso{\code{\link{summary.mat}}, \code{\link{bootstrap.mat}} for boostrap
  resampling of MAT models, \code{\link{predict.mat}} for making
  predictions from MAT models, \code{\link{fitted.mat}} and
  \code{\link{resid.mat}} for extraction of fitted values and residuals
  from MAT models respectively. \code{\link{plot.mat}} provides a
  \code{\link{plot.lm}}-like plotting tool for MAT models.}
\examples{
## Imbrie and Kipp Sea Surface Temperature
data(ImbrieKipp)
data(SumSST)
data(V12.122)

## merge training set and core samples
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- dat[[1]] / 100
ImbrieKippCore <- dat[[2]] / 100

## fit the MAT model using the squared chord distance measure
ik.mat <- mat(ImbrieKipp, SumSST, method = "chord")
ik.mat

## model summary
summary(ik.mat)

## fitted values
fitted(ik.mat)

## model residuals
resid(ik.mat)

## draw summary plots of the model
par(mfrow = c(2,2))
plot(ik.mat)
par(mfrow = c(1,1))

## reconstruct for the RLGH core data
coreV12.mat <- predict(ik.mat, V12.122, k = 3)
coreV12.mat
summary(coreV12.mat)

## draw the reconstruction
reconPlot(coreV12.mat, use.labels = TRUE, display.error = "bars",
          xlab = "Depth", ylab = "SumSST")

## fit the MAT model using the squared chord distance measure
## and restrict the number of analogues we fit models for to 1:20
ik.mat2 <- mat(ImbrieKipp, SumSST, method = "chord", kmax = 20)
ik.mat2

## model summary
summary(ik.mat2)

## fitted values
fitted(ik.mat2)

## model residuals
resid(ik.mat2)

## draw summary plots of the model
par(mfrow = c(2,2))
plot(ik.mat2)
par(mfrow = c(1,1))

## reconstruct for the V12.122 core data
coreV12.mat2 <- predict(ik.mat, V12.122, k = 3)
coreV12.mat2
summary(coreV12.mat2)

## draw the reconstruction
reconPlot(coreV12.mat2, use.labels = TRUE, display.error = "bars",
          xlab = "Depth", ylab = "SumSST")

}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
\keyword{methods}
