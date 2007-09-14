\name{roc}
\alias{roc}
\alias{roc.default}
\alias{roc.mat}
\alias{roc.analog}
\alias{print.roc}
\concept{ROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ROC curve analysis}
\description{
  Fits Receiver Operator Characteristic (ROC) curves to training set
  data. Used to determine the critical value of a dissimilarity
  coefficient that best descriminate between assemblage-types in
  palaeoecological data sets, whilst minimising the false positive
  error rate (FPF).
}
\usage{
roc(object, groups, ...)

\method{roc}{mat}(object, groups, ...)

\method{roc}{analog}(object, groups, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object.}
  \item{groups}{numeric; a vector of group memberships, one entry per
    sample in the training set data.}
  \item{\dots}{arguments passed to/from other methods.}
}
\details{
  A ROC curve is generated from the within-group and between-group
  dissimilarities.

  Within-group dissimilarities are the cells in the lower-triangle
  dissimilarity matrix representing the pairwise dissimilarities for
  samples amples in the same group, over all groups.

  The between-group dissimilarities are the cells in the lower-triangle
  dissimilarity matrix between samples of in a group and all samples not
  in that group, over all groups.
}
\value{
  A list, with the following components:
  
  \item{TPF}{The true positive fraction.}
  \item{FPE}{The false positive error}
  \item{roc.points}{The unique dissimilarities at which the ROC curve
    was evaluated}
  \item{roc.values}{The difference between TPF and FPE at each evaluated
    point of the ROC curve.}
  \item{optimal}{The optimal dissimilarity value, asessed where
    \code{roc.values} is maximal.}
  \item{wilcox}{an object of class \code{"htest"}, the result of a call
    to \code{\link{wilcox.test}}. Contains the results of a Wilcoxon
    Rank Sum and Signed Rank test applied to the within and between
    group dissimilarities.}
  \item{AUC}{The area under the ROC curve.}
  \item{n.within}{numeric; the number of within group dissimilarities.}
  \item{n.without}{numeric; the number of outside of group
    dissimilarities.}
  \item{groups}{numeric; the group membership}
  \item{dissims}{numeric; vector of observed dissimilarities.}
  \item{method}{character; the dissimilarity coefficient used. Taken
    from \code{"object"}.}
  \item{call}{the matched call.}
}
\references{
  
  Brown, C.D., and Davis, H.T. (2006) Receiver operating characteristics
  curves and related decision measures: A tutorial. \emph{Chemometrics
    and Intelligent Laboratory Systems} \bold{80}, 24--38.
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.

  Henderson, A.R. (1993) Assessing test accuracy and its clinical
  consequences: a primer for receiver operating characteristic curve
  analysis. \emph{Annals of Clinical Biochemistry} \strong{30},
  834--846.

}
\author{Gavin L. Simpson, based on code from Thomas Lumley to optimise
  the calculation of the ROC curve.}
\seealso{\code{\link{mat}} for fitting of MAT models.
  \code{\link{bootstrap.mat}} and \code{\link{mcarlo}} for alternative
  methods for selecting critical values of analogue-ness for
  dissimilarity coefficients.}
\examples{
## continue the example from join()
example(join)

## fit a analogue matching (AM) model using the squared chord distance
## measure - need to keep the training set dissimilarities
swap.ana <- analog(swapdiat, rlgh, method = "SQchord",
                   keep.train = TRUE)

## fit the ROC curve to the SWAP diatom data using the AM results
## Generate a grouping for the SWAP lakes
clust <- hclust(as.dist(swap.ana$train), method = "ward")
grps <- cutree(clust, 12)

## fit the ROC curve
swap.roc <- roc(swap.ana, groups = grps)
swap.roc

}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}
\keyword{methods}
