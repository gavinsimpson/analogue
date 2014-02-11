\name{timetrack}
\alias{timetrack}
\alias{print.timetrack}
\alias{plot.timetrack}
\alias{points.timetrack}
\alias{fitted.timetrack}
\alias{scores.timetrack}
\alias{predict.timetrack}

\title{Timetracks of change in species composition}

\description{
  Project passive (e.g. sediment core) samples into an ordination of a
  set of training samples.
}

\usage{
timetrack(X, passive, env, method = c("cca", "rda"),
          transform = "none", formula, scaling = 3,
          rank = "full", \dots)

\method{fitted}{timetrack}(object, which = c("passive", "ordination"),
       model = NULL, choices = 1:2, \dots)

\method{predict}{timetrack}(object, newdata, \dots)

\method{scores}{timetrack}(x, which = c("ordination", "passive"),
       scaling = x$scaling, choices = 1:2, display = "sites", \dots)

\method{plot}{timetrack}(x, choices = 1:2, display = c("wa", "lc"),
     order, type = c("p", "n"), ptype = c("l", "p", "o", "b", "n"),
     pch = c(1,2), col = c("black","red"), lty = "solid", lwd = 1,
     xlim = NULL, ylim = NULL, \dots)

\method{points}{timetrack}(x, choices = 1:2, which = c("passive", "ordination"),
      display = c("wa","lc"), order, \dots)
}

\arguments{
  \item{X}{matrix-like object containing the training set or reference
    samples.}
  \item{passive}{matrix-like object containing the samples to be
    projected into the ordination of \code{X}. Usually a set of sediment
    core samples.}
  \item{env}{optional data frame of environmental or constraining
    variables. If provided, a constrained ordination of \code{X} is
    performed. If \code{formula} is supplied variables named in
    \code{formula} are looked up with \code{env}.}
  \item{method}{character, resolving to an ordination function available
    in \pkg{vegan}. Currently only \code{"cca"}, the default, and
    \code{"rda"} are supported.}
  \item{transform}{character; the name of the transformation to apply to
    both \code{X} and \code{passive}. The transformations are performed
    using \code{tran} and valid options are given by that function's
    \code{method} argument.}
  \item{formula}{a one-sided model formula; if provided, it defines the
    right hand side of the model formula for the ordination function and
    is supplied as argument \code{formula} to the ordination
    function. E.g.~\code{formula = ~ var1 + var2}. If supplied then
    \code{env} must also be supplied}
  \item{scaling}{numeric; the ordination scaling to apply. Useful
    options are likely to be \code{1} or \code{3} where the focus is on
    the samples.}
  \item{rank}{character; see argument of same name in function
    \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.}
  \item{object, x}{an object of class \code{"timetrack"}.}
  \item{which}{character; which fitted values should be returned?}
  \item{model}{character; which ordination component should be used for
    the fitted values; the constrained or unconstrained part? See
    \code{\link{fitted.cca}} for details, but essentially, one of
    \code{"CCA"} for the constrained part and \code{"CA"} for the
    unconstrained part. If \code{NULL}, the default, \code{"CA"} is used
    unless the underlying ordination was constrained, in which case
    \code{"CCA"} is used.}
  \item{choices}{numeric; the length-2 vector of ordination axes to
    plot.}
  \item{newdata}{a data frame of new observations for which locations in
    the plot (or a timetrack) are required. This need not have exactly
    the same set of species as the fitted ordination as internally only
    those species in \code{newdata} that were included in the data used
    for the ordination will be retained. In addition, if a
    transformation was applied to the species data used to fit the
    ordination, the same transformation will be automatically applied to
    \code{newdata} using \code{\link{tran}}.}
  \item{display}{character; which type of sites scores to display? See
    \code{\link{scores.cca}} for details.}
  \item{order}{numeric; vector of indices to use to reorder the passive
    samples. Useful to get passive samples into temporal order for
    plotting with a line.}
  \item{type}{character; the type of plotting required for the training
    set samples. Options are \code{"p"} for points or \code{"n"} to not
    draw training set samples.}
  \item{ptype}{character; controls how the time track should be
    drawn. Default is draw the passive samples connected by a line in
    the order in which they appear in the data. With \code{ptype = "p"}
    no line is drawn. The other types have their usual meaning from
    \code{\link{plot.default}}.}
  \item{pch}{The length-2 vector of plotting characters. The first
    element is used for the ordination samples, the second for the
    passive samples.}
  \item{col}{The length-2 vector of plotting colours. The first
    element is used for the ordination samples, the second for the
    passive samples.}
  \item{lty, lwd}{graphical parameters for the plotted time track for
    \code{ptype != "p"}.}
  \item{xlim, ylim}{user specified axis limits for the plot.}
  \item{\dots}{arguments passed to other methods.
    \code{timetrack} passes arguments on to \code{tran} and the
    ordination function given in \code{method}. \code{fitted} passes
    arguments on to other \code{fitted} methods as
    appropriate. \code{plot} passes arguments on to the underlying
    plotting functions. \code{predict} passes arguments on to
    \code{\link{tran}} for use in applyign the transformation.}
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
  device, whilst the \code{fitted} and \code{scores} methods return the
  matrix of fitted locations on the set of ordination axes.

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
  \item{X}{The training data.}
  \item{transform}{The transformation applied, if any.}
}

\author{
  Gavin L. Simpson
}

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
plot(mod, ptype = "b", col = c("forestgreen", "orange"), lwd = 2)

## Other options (reorder the time track)
ord <- rev(seq_len(nrow(rlgh)))
plot(mod, choices = 2:3, order = ord, ptype = "b",
     col = c("forestgreen", "orange"), lwd = 2)

## illustrating use of the formula
data(swappH)
mod2 <- timetrack(swapdiat, rlgh, env = data.frame(pH = swappH),
                  transform = "hellinger", method = "rda",
                  formula = ~ pH)
mod2
plot(mod2)

## scores and fitted methods
head(fitted(mod, type = "passive"))
head(scores(mod, type = "passive"))

## predict locations in timetrack for new observations
take <- rlgh[1:50, ]
take <- take[ , colSums(take) > 0]
mod3 <- predict(mod, newdata = take)
class(mod3) ## returns a timetrack object
take <- rlgh[-(1:50), ]
take <- take[ , colSums(take) > 0]
mod4 <- predict(mod, newdata = take)

## build a plot up from base parts
plot(mod, type = "n", ptype = "n")
points(mod, which = "ordination", col = "grey", pch = 19, cex = 0.7)
points(mod3, which = "passive", col = "red")
points(mod4, which = "passive", col = "blue")
}

\keyword{methods}
\keyword{hplot}
