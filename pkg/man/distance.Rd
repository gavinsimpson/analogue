\name{distance}
\alias{distance}
\alias{distance.default}
\alias{distance.join}
\concept{dissimilarity}
\concept{dissimilarity coefficient}
\concept{similarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flexibly calculate dissimilarity or distance measures }
\description{
  Flexibly calculates distance or dissimilarity measures between a
  training set \code{x} and a fossil or test set \code{y}. If
  \code{y} is not supplied then the pairwise dissimilarities between
  samples in the training set, \code{x}, are calculated.
}
\usage{

distance(x, \dots)

\method{distance}{default}(x, y, method = c("euclidean", "SQeuclidean",
         "chord", "SQchord", "bray", "chi.square",
         "SQchi.square", "information", "chi.distance",
         "manhattan", "kendall", "gower", "alt.gower",
         "mixed"),
         fast = TRUE,
         weights = NULL, R = NULL, \dots)

\method{distance}{join}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data frame or matrix containing the training set samples, or
    and object of class \code{\link{join}}.}
  \item{y}{data frame or matrix containing the fossil or test set
    samples.}
  \item{method}{character; which choice of dissimilarity coefficient to
    use. One of the listed options. See Details below.}
  \item{fast}{logical; should fast versions of the dissimilarities be
    calculated? See details below.}
  \item{weights}{numeric; vector of weights for each descriptor.}
  \item{R}{numeric; vector of ranges for each descriptor.}
  \item{\dots}{arguments passed to other methods}
}
\details{
  A range of dissimilarity coefficients can be used to calculate
  dissimilarity between samples. The following are currently available:

  \tabular{ll}{
    \code{euclidean}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2}}{d[jk] = sqrt(sum
      (x[ij]-x[ik])^2)}
    \cr
    \code{SQeuclidean}
    \tab \eqn{d_{jk} = \sum_i (x_{ij}-x_{ik})^2}{d[jk] = sum
      (x[ij]-x[ik])^2}
    \cr
    \code{chord}
    \tab \eqn{d_{jk} = \sqrt{\sum_i
	(\sqrt{x_{ij}}-\sqrt{x_{ik}})^2}}{d[jk] = sqrt(sum((sqrt(x) -
      sqrt(y))^2))}
    \cr
    \code{SQchord}
    \tab \eqn{d_{jk} = \sum_i (\sqrt{x_{ij}}-\sqrt{x_{ik}})^2}{d[jk] =
      sum((sqrt(x) - sqrt(y))^2)}
    \cr
    \code{bray}
    \tab \eqn{d_{jk} = \frac{\sum_i |x_{ij} - x_{ik}|}{\sum_i (x_{ij} +
	x_{ik})}}{d[jk] = sum(abs(x - y)) / sum(x + y)}
    \cr
    \code{chi.square}
    \tab \eqn{d_{jk} = \sqrt{\sum_i \frac{(x_{ij} - x_{ik})^2}{x_{ij} +
	x_{ik}}}}{d[jk] = sqrt(sum(((x - y)^2) / (x + y)))}
    \cr
    \code{SQchi.square}
    \tab \eqn{d_{jk} = \sum_i \frac{(x_{ij} - x_{ik})^2}{x_{ij} +
	x_{ik}}}{d[jk] = sum(((x - y)^2) / (x + y))}
    \cr
    \code{information}
    \tab \eqn{d_{jk} = \sum_i (p_{ij}log(\frac{2p_{ij}}{p_{ij} + p_{ik}})
      + p_{ik}log(\frac{2p_{ik}}{p_{ij} + p_{ik}}))}{d[jk] = sum((x[ij] *
      log((2 * x[ij]) / (x[ij] + x[ik]))) + (x[ik] * log((2 * x[ik]) /
      (x[ij] + x[ik]))))}
    \cr
    \code{chi.distance}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2 / (x_{i+} /
	x_{++})}}{d[jk] = sqrt(sum((x[ij] - x[ik])^2 / (x[i+] / x[++])))}
    \cr
    \code{manhattan}
    \tab \eqn{d_{jk} = \sum_i (|x_{ij}-x_{ik}|)}{d[jk] = sum
      (|x[ij]-x[ik]|)}
    \cr
    \code{kendall}
    \tab \eqn{d_{jk} = \sum_i MAX_i - minimum(x_{ij}, x_{ik})}{d[jk] = sum
      (MAX[i] - min(x[ij]-x[ik]))}
    \cr
    \code{gower}
    \tab \eqn{d_{jk} = \sum_i\frac{|p_{ij} -
	  p_{ik}|}{R_i}}{d[jk] = sum(abs(x[ij] - x[ik]) / R[i])}
    \cr
    \code{alt.gower}
    \tab \eqn{d_{jk} = \sqrt{2\sum_i\frac{|p_{ij} -
	  p_{ik}|}{R_i}}}{d[jk] = sqrt(2 * sum(abs(x[ij] - x[ik]) / R[i]))}
    \cr
    \tab where \eqn{R_i}{R[i]} is the range of proportions for
    descriptor (variable) \eqn{i}
    \cr
    \code{mixed}
    \tab \eqn{d_{jk} = \frac{\sum_{i=1}^p w_{i}s_{jki}}{\sum_{i=1}^p
	w_{i}}}{d[jk] = sum(w[i] * s[jki]) / sum(w[i])}
    \cr
    \tab where \eqn{w_i}{w[i]} is the weight for descriptor \eqn{i} and
    \eqn{s_{jki}}{s[jki]} is the similarity \cr
    \tab between samples \eqn{j} and \eqn{k} for descriptor (variable)
    \eqn{i}.
  }

  Argument \code{fast} determines whether fast C versions of some of the
  dissimilarity coefficients are used. The fast versions make use of
  \code{\link[stats]{dist}} for \code{method}s \code{"euclidean"},
  \code{"SQeuclidean"}, \code{"chord"}, \code{"SQchord"}, and
  \code{\link[vegan]{vegdist}} for \code{method} == \code{"bray"}. These
  fast versions are used only when \code{x} is supplied, not when
  \code{y} is also supplied. Future versions of \code{distance} will
  include fast C versions of all the dissimilary coefficients and for
  cases where \code{y} is supplied.
}
\value{
  A matrix of dissimilarities where columns are the samples in
  \code{y} and the rows the samples in \code{x}. If \code{y} is
  not provided then a square, symmetric matrix of pairwise sample
  dissimilarities for the training set \code{x} is returned.

  The dissimilarity coefficient used (\code{method}) is returned as
  attribute \code{"method"}. Attribute \code{"type"} indicates whether
  the object was computed on a single data matrix (\code{"symmetric"})
  or across two matrices (i.e. the dissimilarties between the rows of
  two matrices; \code{"asymmetric"}.
}
\note{
  The dissimilarities are calculated in native R code. As such, other
  implementations (see See Also below) will be quicker. This is done for
  one main reason - it is hoped to allow a user defined function to be
  supplied as argument \code{"method"} to allow for user-extension of
  the available coefficients.

  The other advantage of \code{distance} over other implementations, is
  the simplicity of calculating only the required pairwise sample
  dissimilarities between each fossil sample (\code{y}) and each
  training set sample (\code{x}). To do this in other implementations,
  you would need to merge the two sets of samples, calculate the full
  dissimilarity matrix and then subset it to achieve similar results.
}
\section{warning}{
  For \code{method = "mixed"} it is essential that a factor in \code{x}
  and \code{y} have the same levels in the two data frames. Previous
  versions of analogue would work even if this was not the case, which
  will have generated incorrect dissimilarities for \code{method =
  "mixed"} for cases where factors for a given species had different
  levels in \code{x} to \code{y}. 

  \code{distance} now checks for matching levels for each species
  (column) recorded as a factor. If the factor for any individual
  species has different levels in \code{x} and \code{y}, an error will
  be issued.
}
\references{

  Faith, D.P., Minchin, P.R. and Belbin, L. (1987) Compositional
  dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} \strong{69}, 57--68.
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.

  Kendall, D.G. (1970) A mathematical approach to
  seriation. \emph{Philosophical Transactions of the Royal Society of
    London - Series B} \strong{269}, 125--135.

  Legendre, P. and Legendre, L. (1998) \emph{Numerical Ecology}, 2nd
  English Edition. Elsevier Science BV, The Netherlands.
  
  Overpeck, J.T., Webb III, T. and Prentice I.C. (1985) Quantitative
  interpretation of fossil pollen spectra: dissimilarity coefficients and
  the method of modern analogues. \emph{Quaternary Research} \strong{23},
  87--108. 
  
  Prentice, I.C. (1980) Multidimensional scaling as a research tool in
  Quaternary palynology: a review of theory and methods. \emph{Review of
    Palaeobiology and Palynology} \strong{31}, 71--104.
 
}
\author{Gavin L. Simpson }
\seealso{\code{\link[vegan]{vegdist}} in package \pkg{vegan},
  \code{\link[cluster]{daisy}} in package \pkg{cluster}, and
  \code{\link[stats]{dist}} provide comparable functionality for the
  case of missing \code{y} and are implemented in compiled code, so
  will be faster.}
\examples{
## simple example using dummy data
train <- data.frame(matrix(abs(runif(200)), ncol = 10))
rownames(train) <- LETTERS[1:20]
colnames(train) <- as.character(1:10)
fossil <- data.frame(matrix(abs(runif(100)), ncol = 10))
colnames(fossil) <- as.character(1:10)
rownames(fossil) <- letters[1:10]

## calculate distances/dissimilarities between train and fossil
## samples
test <- distance(train, fossil)

## using a different coefficient, chi-square distance
test <- distance(train, fossil, method = "chi.distance")

## calculate pairwise distances/dissimilarities for training
## set samples
test2 <- distance(train)

## Using distance on an object of class join
dists <- distance(join(train, fossil))
str(dists)

## calculate Gower's general coefficient for mixed data
## first, make a couple of variables factors
fossil[,4] <- factor(sample(rep(1:4, length = 10), 10))
train[,4] <- factor(sample(rep(1:4, length = 20), 20))
## now fit the mixed coefficient
test3 <- distance(train, fossil, "mixed")

## Example from page 260 of Legendre & Legendre (1998)
x1 <- t(c(2,2,NA,2,2,4,2,6))
x2 <- t(c(1,3,3,1,2,2,2,5))
Rj <- c(1,4,2,4,1,3,2,5) # supplied ranges

distance(x1, x2, method = "mixed", R = Rj)

## note this gives 1 - 0.66 (not 0.66 as the answer in
## Legendre & Legendre) as this is expressed as a
## distance whereas Legendre & Legendre describe the
## coefficient as similarity coefficient
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{methods}
