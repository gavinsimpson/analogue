\name{analog}
\alias{analog}
\alias{analog.default}
\alias{analog.distance}
\alias{print.analog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analogue matching}
\description{
  Analogue matching is a  more general implementation of the modern
  analogue methodology than MAT, where we are only interested in
  identifying sufficiently similar samples from a modern training as
  being suitable modern analogues for one or more fossil samples.
}
\usage{
analog(x, \dots)

\method{analog}{default}(x, y,
       method = c("euclidean", "SQeuclidean", "chord", "SQchord",
                  "bray", "chi.square", "SQchi.square",
                  "information", "chi.distance", "manhattan",
                  "kendall", "gower", "alt.gower", "mixed"),
       keep.train = TRUE, \dots)

\method{analog}{distance}(x, train = NULL, keep.train = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{data frames with same columns. \code{x} is training data
    and \code{y}, the test data.}
  \item{method}{character string naming the dissimilarity methods to be
    used. See Details below.}
  \item{keep.train}{logical; should the dissimilarity matrix for the
    training set be stored?}
  \item{train}{a pre-computed dissimilarity matrix for the training set
    samples. Objects of classes \code{"dist"}, \code{"vegdist"}, and
    \code{"distance"} are currently accepted.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  \code{analog} implements analogue matching \emph{sensu} Flower et
  al (1997) and Simpson et al (2005), where the aim is to identify
  suitable close analogues of fossil samples from a modern training
  set. These results are generally used within ecological restoration,
  but the identification of close modern analogues for fossil samples is
  also used as a technique for assessing transfer function
  reconstructions.

  \code{analog} is a simple and very general function that generates a
  pairwise dissimilarity matrix for the modern training set, and a second
  matrix containing the pairwise dissimilarities between each fossil
  sample and each sample in the training set. These results can then be
  assessed using other functions and to extract the close modern
  analogues using function \code{\link{cma}}. See the See Also section
  below.

  Analysis of the pairwise dissimilarity matrix for the modern training
  set can be used to decide on a suitable dissimilarity threshold
  for defining close modern analogues. By default this matrix is
  returned as part of the output from the \code{analog} function.
}
\value{
  A list of class \code{"analog"} with the following components:
  \item{analogs}{matrix of pairwise dissimilarities between each fossil
    sample (\code{y}) and each sample in the modern training set
    (\code{x}).}
  \item{train}{if argument \code{keep.train} is \code{TRUE} then a
    pairwise dissimilarity matrix for the modern training set.}
  \item{call}{the matched function call.}
  \item{method}{character; the dissimilarity coefficient used.}
}
\references{

  Flower, R.J., Juggins, S. and Battarbee, R.W. (1997) Matching diatom
  assemblages in lake sediment cores and modern surface sediment
  samples: the implications for lake conservation and restoration with
  special reference to acidified systems. \emph{Hydrobiologia}
  \strong{344}; 27--40.

  Simpson, G.L., Shilland, E.M., Winterbottom, J. M. and Keay, J. (2005)
  Defining reference conditions for acidified waters using a modern
  analogue approach. \emph{Environmental Pollution} \strong{137};
  119--133.
  
}
\author{Gavin L. Simpson}
\seealso{
  \code{\link{distance}} for the function that calculates the
  dissimilarity matrices.
  \code{\link{cma}} for extraction of close modern analogues.
  \code{\link{dissimilarities}} and \code{\link{plot.dissimilarities}}
  for analysis of distribution of pairwise dissimilarity matrix for
  modern training set.
}
\examples{
## Imbrie and Kipp example
## load the example data
data(ImbrieKipp)
data(SumSST)
data(V12.122)

## merge training and test set on columns
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- dat[[1]] / 100
V12.122 <- dat[[2]] / 100

## Imbrie and Kipp foraminfera sea-surface temperature

## analog matching between SWAP and RLGH core
ik.analog <- analog(ImbrieKipp, V12.122, method = "chord")
ik.analog
summary(ik.analog)

## Can take pre-computed dissimilarity objects
d1 <- distance(ImbrieKipp, V12.122)
d2 <- distance(ImbrieKipp)
ik <- analog(d1, d2, keep.train = TRUE)
ik

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
