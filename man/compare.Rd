\name{compare}
\alias{compare}
\alias{compare.default}

\title{Compare proxies across two data sets}
\description{
  \code{compare()} compares a proxy dataset with a training set or other
  data set that is considered to be the master. A range of metrics is
  returned, currently for samples only.
}
\usage{
compare(x, \dots)

\method{compare}{default}(x, y, env,
        by = c("sites", "species"),
        ordination = "rda",
        method = "chord",
        transform = NULL,
        n2limit = 5L,
        ...)
}

\arguments{
  \item{x}{data frame; training set samples to compare against}
  \item{y}{data frame; passive or core samples}
  \item{env}{numeric vector of environmental or contraint data for
    residual length ordination. Ignored if \code{by = "species"}.}
  \item{by}{character; compare data sets by sites or species
    (proxies).}
  \item{ordination}{character; which constrained ordination method to
    use}
  \item{method}{character; which dissimilarity method to use. See
    \code{distance}.}
  \item{transform}{character: should a transformation be applied to the
    data. Ignored.}
  \item{n2limit}{integer; the criterion for indicating species with
    potentially poorly estimated optima. The default value of \code{5L}
    is one suggested by R. Telford.}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  ToDo
}
\value{
  If \code{by = "species"} a data frame of diagnostics for each species
  (proxy) in \code{y} relative to \code{x}. If \code{by = "sites"}, the
  diagnostics are for each sample (row) in \code{y}. Depending on the
  value of \code{by} some of the following columns will be returned

  \item{\code{sumMissing}}{numeric; abundance sum for species missing
    from the training set \code{x}.}
  \item{\code{sumPoorOpt}}{numeric; abundance sum for species with
    potentially poorly estimated optima.}
  \item{\code{closestSamp} }{numeric; minimum dissimilarity to a sample
    in the training data \code{x}.}
  \item{\code{residLen}}{numeric; the squared residual length for each
    sample in \code{y}. A measure of how well the sample fits into the
    species-environment relationship of a constrained ordination. See
    \code{\link{residLen}} for details. Not given if \code{env} is not
    provided.}
  \item{\code{inTrain}}{logical; simple indicator of whether a species
    in \code{y} is present in the training data \code{x}.}
  \item{\code{n2}}{numeric; Hill's N2 for each species in \code{y}.}
  \item{\code{n2Train}}{numeric; as for \code{n2} but computed from the
    training data \code{x}.}
  \item{\code{max}}{numeric; the maximum abundance of each species
    computed using \code{y}.}
  \item{\code{maxTrain}}{numeric; as for \code{max} but computed using
    the training data \code{x}.}
}

\author{Gavin L. Simpson}

\examples{
data(ImbrieKipp, V12.122, SumSST)
compare(ImbrieKipp, V12.122, env = SumSST, ordination = "rda",
        method = "chord")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ utility }
\keyword{ multivariate }
