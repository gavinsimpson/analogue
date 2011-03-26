\name{roc}
\alias{roc}
\alias{roc.default}
\alias{roc.mat}
\alias{roc.analog}
\alias{print.roc}
\alias{summary.roc}
\alias{print.summary.roc}
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
roc(object, groups, k = 1, ...)

\method{roc}{default}(object, groups, k = 1, thin = FALSE,
    max.len = 10000, ...)

\method{roc}{mat}(object, groups, k = 1, ...)

\method{roc}{analog}(object, groups, k = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object.}
  \item{groups}{a vector of group memberships, one entry per
    sample in the training set data. Can be a factor, and will be
    coerced to one if supplied vecvtor is not a factor.}
  \item{k}{numeric; the \code{k} closest analogues to use to calculate
    ROC curves.}
  \item{thin}{logical; should the points on the ROC curve be thinned?
    See Details, below.}
  \item{max.len}{numeric; length of analolgue and non-analogue
    vectors. Used as limit to thin points on ROC curve to.}
  \item{\dots}{arguments passed to/from other methods.}
}
\details{
  A ROC curve is generated from the within-group and between-group
  dissimilarities.

  For each level of the grouping vector (\code{groups}) the
  dissimilarity between each group member and it's k closest analogues
  within that group are compared with the k closest dissimilarities
  between the non-group member and group member samples.

  If one is able to discriminate between members of different group on
  the basis of assemblage dissimilarity, then the dissimilarities
  between samples within a group will be small compared to the
  dissimilarities between group members and non group members.

  \code{thin} is useful for large problems, where the number of analogue
  and non-analogue distances can conceivably be large and thus overflow
  the largest number R can work with. This option is also useful to
  speed up computations for large problems. If \code{thin == TRUE}, then
  the larger of the analogue or non-analogue distances is thinned to a
  maximum length of \code{max.len}. The smaller set of distances is
  scaled proportionally. In thinning, we approximate the distribution of
  distances by taking \code{max.len} (or a fraction of \code{max.len}
  for the smaller set of distances) equally-spaced probability
  quantiles of the distribution as a new set of distances. 
}
\value{
  A list with two components; i, \code{statistics}, a summary of ROC
  statistics for each level of \code{groups} and a combined ROC
  analysis, and ii, \code{roc}, a list of ROC objects, one per level of
  \code{groups}. For the latter, each ROC object is a list, with the
  following components: 
  
  \item{TPF}{The true positive fraction.}
  \item{FPE}{The false positive error}
  %\item{roc.values}{The difference between TPF and FPE at each evaluated
  %  point of the ROC curve.}
  \item{optimal}{The optimal dissimilarity value, asessed where the
  slope of the ROC curve is maximal.}
  \item{AUC}{The area under the ROC curve.}
  \item{se.fit}{Standard error of the AUC estimate.}
  \item{n.in}{numeric; the number of samples within the current group.}
  \item{n.out}{numeric; the number of samples not in the current group.}
  \item{p.value}{The p-value of a Wilcoxon rank sum test on the two sets
    of dissimilarities. This is also known as a Mann-Whitney test.}
  \item{roc.points}{The unique dissimilarities at which the ROC curve
    was evaluated}
  \item{max.roc}{numeric; the position along the ROC curve at which the
    slope of the ROC curve is maximal. This is the index of this point
    on the curve.}
  \item{prior}{numeric, length 2. Vector of observed prior probabilities
    of true analogue and true non-analogues in the group.}
  \item{analogue}{a list with components \code{yes} and \code{no}
    containing the dissimilarities for true analogue and true
    non-analogues in the group.}
  %\item{groups}{numeric; the group membership}
  %\item{dissims}{numeric; vector of observed dissimilarities.}
  %\item{method}{character; the dissimilarity coefficient used. Taken
  %  from \code{"object"}.}
  %\item{call}{the matched call.}
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
## load the example data
data(swapdiat, swappH, rlgh)

## merge training and test set on columns
dat <- join(swapdiat, rlgh, verbose = TRUE)

## extract the merged data sets and convert to proportions
swapdiat <- dat[[1]] / 100
rlgh <- dat[[2]] / 100

## fit an analogue matching (AM) model using the squared chord distance
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

## draw the ROC curve
plot(swap.roc, 1)

## draw the four default diagnostic plots
layout(matrix(1:4, ncol = 2))
plot(swap.roc)
layout(1)
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}
\keyword{methods}
