\name{evenSample}
\alias{evenSample}

\title{Number of samples per gradient segments}

\description{
  The number of samples in sections along the gradient is a useful
  diagnostic as to the quality of reconstructions at gradient values
  within those sections.
}

\usage{
evenSample(grad, n = 10)
}

\arguments{
  \item{grad}{numeric; vector of gradient values}
  \item{n}{number of segments to partition the gradient into}
}

\details{
  The sampling design of a training set, i.e. the number of samples taken
  at points along the gradient, can influence the uncertainty in the
  transfer function predictions at those values of the gradient. Poorly
  sampled sections of the gradient may have far larger RMSEP than the
  overall model RMSEP.
}

\value{
  Numeric vector of length \code{n} containing the numbers of samples
  per gradient segment.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{Gavin L. Simpson}

\examples{
data(SumSST)
evenSample(SumSST) ## not an even sample...
}

\keyword{ utilities }
