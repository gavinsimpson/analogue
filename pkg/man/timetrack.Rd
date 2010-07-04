\name{timetrack}
\alias{timetrack}
\alias{print.timetrack}
\alias{plot.timetrack}
\alias{fitted.timetrack}

\title{Timetracks of change in species composition}

\description{
  Project passive (e.g. sediment core) samples into an ordination of a
  set of training samples.
}

\usage{
timetrack(X, passive, env, method = c("cca", "rda"),
          transform = "none", formula, scaling = 3,
          rank = "full", model = c("CCA", "CA"), \dots)

\method{fitted}{timetrack}(object, type = c("passive", "ordination"),
       model = NULL, \dots)

\method{plot}{timetrack}(x, choices = 1:2, pch = c(1,2),
     col = c("black","red"), \dots)
}

\arguments{
  \item{X}{matrix-like object containing the training set or reference
    samples.}
  \item{passive}{matrix-like object containing the samples to be
    projected into the ordination of \code{X}. Usually a set of sediment
    core samples.}
  \item{env}{optional vector or matrix of environmental or constraining
    variables. If provided, a constrained ordination of \code{X} is
    performed.}
  \item{method}{character, resolving to an ordination function available
    in \pkg{vegan}. Currently only \code{"cca"}, the default, and
    \code{"rda"} are supported.}
  \item{transform}{character; the name of the transformation to apply to
    both \code{X} and \code{passive}. The transformations are performed
    using \code{tran} and valid options are given by that function's
    \code{method} argument.}
  \item{formula}{a model formula; if provided, it defines the model
    formula for the ordination function and is supplied as argument
  \code{formula} to the ordination function.}
  \item{scaling}{numeric; the ordination scaling to apply. Useful
    options are likely to be \code{1} or \code{3} where the focus is on
    the samples.}
  \item{rank}{character; see argument of same name in function
    \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.}
  \item{model}{character; see argument of same name in function
    \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.}
  \item{object, x}{an object of class \code{"timetrack"}.}
  \item{type}{character; which fitted values should be returned?}
  \item{choices}{numeric; the length-2 vector of ordination axes to
    plot.}
  \item{pch}{The length-2 vector of plotting characters. The first
    element is used for the ordination samples, the second for the
    passive samples.}
  \item{col}{The length-2 vector of plotting colours. The first
    element is used for the ordination samples, the second for the
    passive samples.}
  \item{\dots}{arguments passed to other methods.
    \code{timetrack} passes arguments on to \code{tran} and the
    ordination function given in \code{method}. \code{fitted} passes
    arguments on to other \code{fitted} methods as
    appropriate. \code{plot} passes arguments on to the underlying
    plotting functions.}
}

\details{
  The timetrack is a way to visualise changes in species composition
  from sediment core samples within an underlying reference ordination
  or, usually, training set samples. This technique has been most often
  applied in situations where the underlying ordination is a constrained
  ordination and thence the timetrack of sediment core samples within
  the ordination reflects both the change in species composition and the
  indicative changes in the constraining variables.

  The sediment core samples are projected passively into the underlying
  ordination. By projected passively, the locations of the core samples
  are predicted on the basis of the ordination species scores. A common
  set of species (columns) is required to passively place the sediment
  samples into the ordination. To achieve this, the left outer join of
  the species compositions of the training set and passive set
  is determined; the left outer join results in the passive data matrix
  having the same set of species (variables; columns) as the training
  set. Any training set species not in the passive set are added to
  the passive set with abundance 0. Any passive species not in the
  training set are removed from the passive set.
}

\value{
  The \code{plot} method results in a plot on the currently active
  device, whilst the \code{fitted} method returns the matrix of fitted
  locations on the set of ordination axes.

  \code{timetrack} returns an object of class \code{"timetrack"}, a list
  with the following components:
  \item{ordination }{the ordination object, the result of the call to
    the function of the name \code{method}.}
  \item{fitted.values }{the matrix of fitted locations for the passive
    samples on the ordination axes.}
  \item{method }{the ordination function used.}
  \item{formula }{if supplied, the model formula used to define the
    ordination model.}
  \item{scaling }{the ordination scaling applied.}
  \item{rank }{The rank or the number of axes used in the
    approximation. The default is to use all axes (full rank) of the
    \code{"model"}.} 
  \item{model }{Show constrained (\code{"CCA"}) or unconstrained
    (\code{"CA"}) results.}
  \item{labels }{a list of names for the \code{X}, \code{passive}, and
    \code{env} arguments.}
  \item{call }{The matched function call.}
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Gavin L. Simpson
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}} for the
  underlying ordination functions.
}

\examples{
## load the RLGH and SWAP data sets
data(rlgh, swapdiat)

## Fit the timetrack ordination
mod <- timetrack(swapdiat, rlgh, transform = "hellinger",
                 method = "rda")
mod

## Plot the timetrack
plot(mod)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
