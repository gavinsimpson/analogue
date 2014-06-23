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

\method{compare}{default}(x, y, env, ordination = "rda",
        method = "chord", transform = NULL, ...)
}

\arguments{
  \item{x}{data frame; training set samples to compare against}
  \item{y}{data frame; passive or core samples}
  \item{env}{numeric vector of environmental or contraint data for
    residual length ordination}
  \item{ordination}{character; which constrained ordination method to
    use}
  \item{method}{character; which dissimilarity method to use. See
    \code{distance}.}
  \item{transform}{character: should a transformation be applied to the
    data. Ignored.}
  \item{\dots}{arguments passed to other methods}
}
\details{
  ToDo
}
\value{
  Returns a data frame of diagnostics for each sample in \code{y}.
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
