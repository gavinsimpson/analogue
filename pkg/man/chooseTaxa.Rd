\name{chooseTaxa}
\alias{chooseTaxa}
\alias{chooseTaxa.default}
\title{Select taxa (variables) on basis of maximum abundance attained
  and number of occurrences}
\description{
  Select taxa (variables) from an object on the basis of one or both of
  maximum abundance and number of occurrences greater than user-specified
  values. This is a simple utility function to encapsulate this common
  task in filtering palaeoecological data sets.
}
\usage{
chooseTaxa(object, \dots)

\method{chooseTaxa}{default}(object, n.occ = 1, max.abun = 0,
           type = c("AND","OR"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object for which a suitable method exists. The
    default method assumes a matrix-like object such as a data
    frame or a numeric matrix.}
  \item{n.occ}{numeric; number of occurrences representing the lower
    limit for selection. A taxon is included in the returned subset if
    it is present a total of \code{n.occ} times or more. See argument
    \code{type} for a modifier which might exclude the taxon even if it
    would be included on the basis of \code{n.occ}.}
  \item{max.abun}{numeric; maximum abundance representing the lower
    limit for selection. A taxon is included in the returned subset if
    it attains abundance equal to or greater than \code{max.abun} in one
    or more sample. See argument \code{type} for a modifier which might
    exclude the taxon even if it would be included on the basis
    of \code{max.abun}.}
  \item{type}{character; one of \code{"AND"} or \code{"OR"}, controlling
    how the criteria \code{n.occ} and \code{max.abun} are combined to
    generate a subset of the variables in \code{object}.}
  \item{\dots}{arguments passed on to subsequent methods.}
}
%\details{
%}
\value{
  Returns the supplied data frame or matrix with a subset of columns
  (taxa) that meet the criteria chosen.
}
%\references{
%}
\author{Gavin L. Simpson}
%\seealso{
%}
\examples{
data(ImbrieKipp)
IK2 <- chooseTaxa(ImbrieKipp, n.occ = 5)
dim(ImbrieKipp)
dim(IK2)
}
\keyword{methods}
\keyword{multivariate}
