\name{weightedCor}
\alias{weightedCor}
\alias{weightedCor.default}
\alias{print.weightedCor}
\alias{plot.weightedCor}
\title{
  Weighted correlation test of WA reconstruction
}
\description{
  Weighted correlation between WA optima from training set and axis 1
  scores of constrained ordination fitted to fossil data with WA model
  predictions for fossil samples as constraints.
}
\usage{
\method{weightedCor}{default}(x, env, fossil, method = c("rda", "cca"),
            test = TRUE, type = c("simulate", "permute"), sim = 999,
            verbose = TRUE, ...)

\method{plot}{weightedCor}(x,
     type = c("bubble", "null"),
     weighted = TRUE,
     size = 0.25,
     xlab = paste(x$env, "WA Optima"),
     ylab = "Axis 1 Score",
     xlim,
     main = "",
     sub = NULL,
     border = "gray75",
     col = "gray75",
     obscol = "red",
     fg = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{training set covariates, a matrix-like object usually of
    species/proxy data. For the \code{plot} method, an object of class
    \code{"weightedCor"}, the result of a call to \code{weightedCor}.}
  \item{env}{training set response, a vector usually of environmental
    data.}
  \item{fossil}{matrix of fossil/core species/proxy data for which a
    reconstruction is sought.}
  \item{method}{constrained ordination method. One of \code{"rda"} and
    \code{"cca"}. Currently only \code{"rda"} is supported.}
  \item{test}{logical; should the observed correlation be tested?}
  \item{type}{the type of test to apply. One of \code{"simulate"} or
    \code{"permute"}. The latter is currently not implemented. For the
    \code{plot} method, the type of plot to produce.}
  \item{sim}{numeric; number of simulations or permutations to permform
    as part of the test}
  \item{verbose}{logical; should the progress of the test be shown via a
    progress bar?}
  \item{\dots}{arguments passed to other methods. In the case of the
    \code{plot} method, additional graphical parameters can be
    supplied.}
  \item{weighted}{logical; should the null distribution plotted be of
     the weighted or normal correlation.}
  \item{size}{numeric; the size of the largest bubble in inches. See
    \code{\link{symbols}} and argument \code{inches} for details.}
  \item{xlim,xlab,ylab,main,sub}{graphical parameters with their
    usual meaning.}
  \item{border, col}{The border and fill colours for the histogram
    bars.}
  \item{fg}{The colour of the bubbles drawn on the bubble plot.}
  \item{obscol}{The colour of the indicator for the observed
    correlation.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  The \code{plot} method produces a plot on the current
  device. \code{weightedCor()} returns a list with the following
  components:

  \item{wtdCorrel,Correl}{numeric; the observed weighted and standard
    correlation.}
  \item{data}{data frame; containing the training set WA Optima, axis 1
    species scores, and mean abundance for each species.}
  \item{ord}{the fitted constrained ordination.}
  \item{model}{the fitted WA model.}
  \item{method}{the ordination method used.}
  \item{ndist}{the null distribution produced. \code{NULL} if argument
    \code{test} was \code{FALSE}.}
  \item{sim}{numeric; the number of simulations or permutations used to
    test the observed correlations.}
  \item{type}{the type of test performed.}
  \item{env}{the deparsed version of \code{env} argument. Used for
    plotting.}
  \item{call}{the matched function call.}
}
\references{
  Telford R.J. and Birks, H.J.B. (2011) A novel method for assessing the
  statistical significance of quantitative reconstructions inferred from
  biotic assemblages. \emph{Quanternary Science Reviews}
  \strong{30}:1272-1278.
}

\author{
  Gavin L. Simpson
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{wa}} for details on fitting weighted average models.
}
\examples{
data(ImbrieKipp, SumSST, V12.122)

Cor <- weightedCor(ImbrieKipp, env = SumSST,
                   fossil = V12.122, type = "simulate", sim = 49)
Cor

plot(Cor)
plot(Cor, type = "null")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
