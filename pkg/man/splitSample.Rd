\name{splitSample}
\alias{splitSample}

\title{
  Select samples from along an environmental gradient
}
\description{
  Select samples from along an environmental gradient by splitting
  the gradient into discrete chunks and sample within each chunk. This
  allows a test set to be selected which covers the environmental
  gradient of the training set, for example.
}
\usage{
splitSample(env, chunk = 10, take, nchunk,
            fill = c("head", "tail", "random"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{numeric; vector of samples representing the gradient values.}
  \item{chunk}{numeric; number of chunks to split the gradient into.}
  \item{take}{numeric; how many samples to take from the gradient. Can
    not be missing.}
  \item{nchunk}{numeric; number of samples per chunk. Must be a vector
    of length \code{chunk} and \code{sum(chunk)} must equal
    \code{take}. Can be missing (the default), in which case some simple
    heuristics are used to determine the number of samples chosen per
    chunk. See Details.}
  \item{fill}{character; the type of filling of chunks to perform. See
    Details.}
}
\details{
  The gradient is split into \code{chunk} sections and samples are
  selected from each chunk to result in a sample of length
  \code{take}. If \code{take} is divisible by \code{chunk} without
  remainder then there will an equal number of samples selected from
  each chunk. Where \code{chunk} is not a multiple of \code{take} and
  \code{nchunk} is not specified then extra samples have to be allocated
  to some of the chunks to reach the required number of samples
  selected.

  An additional complication is that some chunks of the gradient may
  have fewer than \code{nchunk} samples and therefore more samples need
  to be selected from the remaining chunks until \code{take} samples are
  chosen.

  If \code{nchunk} is supplied, it must be a vector stating exactly how
  many samples to select from each chunk. If \code{chunk} is not
  supplied, then the number of samples per chunk is determined as
  follows:

  \enumerate{
    \item An intial allocation of \code{floor(take / chunk)} is assigned
    to each chunk
    \item If any chunks have fewer samples than this initial allocation,
    these elements of \code{nchunk} are reset to the number of  samples
    in those chunks
    \item Sequentially an extra sample is allocated to each chunk with
    sufficient available samples until \code{take} samples are
    selected.
  }

  Argument \code{fill} controls the order in which the chunks are
  filled. \code{fill = "head"} fills from the low to the high end of the
  gradient, whilst \code{fill = "tail"} fills in the opposite
  direction. Chunks are filled in random order if \code{fill =
  "random"}. In all cases no chunk is filled by more than one extra
  sample until all chunks that can supply one extra sample are
  filled. In the case of \code{fill = "head"} or \code{fill = "tail"}
  this entails moving along the gradient from one end to the other
  allocating an extra sample to available chunks before starting along
  the gradient again. For \code{fill = "random"}, a random order of
  chunks to fill is determined, if an extra sample is allocated to each
  chunk in the random order and \code{take} samples are still not
  selected, filling begins again using the same random ordering. In
  other words, the random order of chunks to fill is chosen only once.
}
\value{
  A numeric vector of indices of selected samples. This vector has
  attribute \code{lengths} which indicates how many samples were
  actually chosen from each chunk.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Gavin L. Simpson
}

\examples{
data(swappH)

## take a test set of 20 samples along the pH gradient
test1 <- splitSample(swappH, chunk = 10, take = 20)
test1
swappH[test1]

## take a larger sample where some chunks don't have many samples
## do random filling
set.seed(3)
test2 <- splitSample(swappH, chunk = 10, take = 70, fill = "random")
test2
swappH[test2]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{utilities}
