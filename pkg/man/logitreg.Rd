\name{logitreg}
\alias{logitreg}
\alias{logitreg.default}
\alias{logitreg.analog}
\alias{print.logitreg}
\alias{summary.logitreg}
\alias{print.summary.logitreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Logistic regression models for assessing analogues/non-analogues}
\description{
  Fits logistic regression models to each level of \code{group} to
  model the probability of two samples being analogues conditional upon
  the dissimilarity between the two samples.
}
\usage{
logitreg(object, groups, k = 1, biasReduced = FALSE, ...)

\method{logitreg}{default}(object, groups, k = 1, ...)

\method{logitreg}{analog}(object, groups, k = 1, ...)

\method{summary}{logitreg}(object, p = 0.9, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{for \code{logitreg}; a full dissimilarity matrix. For
    \code{summary.logitreg} an object of class \code{"logitreg"}, the
    result of a call to \code{logitreg}.}
  \item{groups}{factor (or object that can be coerced to one) containing
  the group membership for each sample in \code{object}.}
  \item{k}{numeric; the \code{k} closest analogues to use in the model
    fitting.}
  \item{biasReduced}{logical; should Firth's method for bias reduced
    logistic regression be used to fit the models? If \code{TRUE}, model
    fits are performed via \code{\link{brglm}}. The default,
    \code{FALSE}, indicates that models will be fitted via the standard
    \code{\link{glm}} function.}
  \item{p}{probability at which to predict the dose needed.}
  \item{\dots}{arguments passed to other methods. These arguments are
    passed on to \code{\link{glm}} or \code{\link{brglm}}. See their
    respective helps pages for details. Note that \code{logitreg} sets
    internally the \code{formula}, \code{data}, and \code{family}
    arguments and hence can not be specified by the user.}
}
\details{
  Fits logistic regression models to each level of \code{group} to
  model the probability of two samples being analogues (i.e. in the same
  group) conditional upon the dissimilarity between the two samples.

  This function can be seen as a way of directly modelling the
  probability that two sites are analogues, conditional upon
  dissimilarity, that can also be done less directly using
  \code{\link{roc}} and \code{\link{bayesF}}.

  Often, the number of true analogues in the training set is small, both
  in absolute terms and as a proportion of comparisons. Logistic
  regression is known to suffer from a small-sample bias. Firth's method
  of bias reduction is a general solution to this problem and is
  implemented in \code{logitreg} through the \pkg{brglm} package of
  Ioannis Kosmidis.
}
\value{
  \code{logitreg} returns an object of class \code{"logitreg"}; a list
  whose components are objects returned by \code{\link{glm}}. See
  \code{\link{glm}} for further details on the returned objects.

  The components of this list take their names from \code{group}.

  For \code{summary.logitreg} an object of class
  \code{"summary.logitreg"}, a data frame with summary statistics of the
  model fits. The components of this data frame are:

  \item{In, Out}{The number of analogue and non-analogue dissimilarities
    analysed in each group,}
  \item{Est.(Dij), Std.Err}{Coefficient and its standard error for
    dissimilarity from the logit model,}
  \item{Z-value, p-value}{Wald statistic and associated p-value for each
    logit model.}
  \item{Dij(p=?), Std.Err(Dij)}{The dissimilarity at which the posterior
    probability of two samples being analogues is equal to \eqn{p}, and
    its standard error. These are computed using
    \code{\link[MASS]{dose.p}}.}
}
\references{
  Firth, D. (1993). Bias reduction of maximum likelihood
  estimates. \emph{Biometrika} \strong{80}, 27-38.
}
\author{Gavin L. Simpson}
\note{
  The function may generate warnings from function
  \code{\link{glm.fit}}. These should be investigated and not simply
  ignored.

  If the message is concerns fitted probabilities being numerically 0 or
  1, then check the fitted values of each of the models. These may well
  be numerically 0 or 1. Heed the warning in \code{\link{glm}} and read
  the reference cited therein which \strong{may} indicate problems with
  the fitted models, such as (quasi-)complete separation.
}
\seealso{\code{\link{roc}}, \code{\link{bayesF}}, \code{\link{glm}}, and
  \code{\link{brglm}}.}
\examples{
## load the example data
data(swapdiat, swappH, rlgh)

## merge training and test set on columns
dat <- join(swapdiat, rlgh, verbose = TRUE)

## extract the merged data sets and convert to proportions
swapdiat <- dat[[1]] / 100
rlgh <- dat[[2]] / 100

## fit an analogue matching (AM) model using the squared chord distance
## measure - need to keep the training set dissimilarities
swap.ana <- analog(swapdiat, rlgh, method = "SQchord",
                   keep.train = TRUE)

## fit the ROC curve to the SWAP diatom data using the AM results
## Generate a grouping for the SWAP lakes
clust <- hclust(as.dist(swap.ana$train), method = "ward")
grps <- cutree(clust, 6)

## fit the logit models to the analog object
swap.lrm <- logitreg(swap.ana, grps)
swap.lrm

## summary statistics
summary(swap.lrm)

## plot the fitted logit curves
plot(swap.lrm, conf.type = "polygon")

## extract fitted posterior probabilities for training samples
## for the individual groups
fit <- fitted(swap.lrm)
head(fit)

## compute posterior probabilities of analogue-ness for the rlgh
## samples. Here we take the dissimilarities between fossil and
## training samples from the `swap.ana` object rather than re-
## compute them
pred <- predict(swap.lrm, newdata = swap.ana$analogs)
head(pred)

## Bias reduction
## fit the logit models to the analog object
swap.brlrm <- logitreg(swap.ana, grps, biasReduced = TRUE)
summary(swap.brlrm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
