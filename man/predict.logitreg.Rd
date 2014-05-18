\name{predict.logitreg}
\alias{predict.logitreg}

\title{Posterior probability of analogue-ness for fossil samples}

\description{
  Predict the posterior probability of analogue-ness between fossil and
  training set samples based on logistic regression fits.
}

\usage{

\method{predict}{logitreg}(object, newdata, group = "all", k = 1, ...)

}

\arguments{
  \item{object}{an object of class \code{"logitreg"}, resulting from a
    call to \code{\link{logitreg}}.}
  \item{newdata}{an object of class \code{"distance"} representing the
    dissimilarity between fossil samples and those from the training
    set.}
  \item{group}{character; for which group(s) are predictions
    sought. \code{"all"}, the default, indicates that predictions for
    all groups and the combined analysis are produced. Currently
    ignored.}
  \item{k}{numeric; the number of close modern analogues to each fossil
    sample to consider. Not currently used and may be removed from the
    method as varying \code{k} should require refitting the logistic
    regression model with that number of close modern analogues
    considered.}
  \item{\dots}{additional arguments passed to other methods.}
}
% \details{

% }

\value{
  A matrix of posterior probabilities of analogue-ness is returned. The
  rows represent the fossil samples and the columns the groupings. There
  are \eqn{g+1} columns where \eqn{g} is the number of groups. The
  \eqn{g+1}th column represents the Combined analysis which is the
  overall posterior probability that a fossil sample is an analogue to a
  training set sample (i.e. to any one of the groups).
}

\author{Gavin L. Simpson}

\seealso{
  See \code{\link{logitreg}} for example usage. \code{\link{cma}} for
  extraction of close modern analogue information and
  \code{\link{analog}} for the underlying computations.
}

% \examples{

% }

\keyword{methods}
