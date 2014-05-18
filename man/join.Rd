\name{join}
\alias{join}
\alias{head.join}
\alias{tail.join}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge species data sets on common columns (species)}
\description{
  Merges any number of species matrices on their common columns to
  create a new data set with number of columns equal to the number of
  unqiue columns across all data frames. Needed for analysis of fossil
  data sets with respect to training set samples. 
}
\usage{
join(\dots, verbose = FALSE, na.replace = TRUE, split = TRUE, value = 0,
     type = c("outer", "left", "inner"))

\method{head}{join}(x, \dots)

\method{tail}{join}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{for \code{join}, data frames containing the data sets to
    be merged. For the \code{\link{head}} and \code{\link{tail}}
    methods, additional arguments to \code{\link{head}} and
    \code{\link{tail}}, in particular \code{"n"} to control the
    number of rows of each \code{join}ed data set to display.} 
  \item{verbose}{logical; if \code{TRUE}, the function prints out the
    dimensions of the data frames in \code{"\dots"}, as well as those of
    the returned, merged data frame.}
  \item{na.replace}{logical; samples where a column in one data frame
    that have no matching column in the other will contain missing
    values (\code{NA}). If \code{na.replace} is \code{TRUE}, these
    missing values are replaced with zeros. This is standard practice in
    ecology and palaeoecology. If you want to replace with another
    value, then set \code{na.replace} to \code{FALSE} and do the
    replacement later.}
  \item{split}{logical; should the merged data sets samples be split
    back into individual data frames, but now with common columns
    (i.e. species)?}
  \item{value}{numeric; value to replace \code{NA} with if
    \code{na.replace} is \code{TRUE}.}
  \item{type}{logical; type of join to perform. \code{"outer"} returns
    the \emph{union} of the variables in data frames to be merged, such
    that the resulting objects have columns for all variables found
    across all the data frames to be merged. \code{"left"} returns the
    left outer (or the left) join, such that the merged data frames
    contain the set of variables found in the first supplied data
    frame. \code{"inner"} returns the inner join, such that the merged
    data frame contain the intersection of the variables in the supplied
    data frames. See Details.}
  \item{x}{an object of class \code{"join"}, usually the result of a
    call to \code{\link{join}}.}
}
\details{
  When merging multiple data frames the set of variables in the merged
  data can be determined via a number of routes. \code{join} provides
  for two (currently) join types; the \emph{outer} join and the
  \emph{left outer} (or simply the \emph{left}) join. Which type of join
  is performed is determined by the argument \code{type}.

  The \emph{outer} join returns the union of the set of variables found
  in the data frames to be merged. This means that the resulting data
  frame(s) contain columns for all the variable observed across all the
  data frames supplied for merging.

  With the \emph{left outer} join the resulting data frame(s) contain
  only the set of variables found in the first data frame provided.

  The \emph{inner} join returns the intersection of the set of variables
  found in the supplied data frames. The resulting data frame(s)
  contains the variables common to all supplied data frames.
}
\value{
  If \code{split = TRUE}, an object of class \code{"join"}, a list of
  data frames, with as many components as the number of data frames
  originally merged.

  Otherwise, an object of class \code{c("join", "data.frame")}, a data
  frame containing the merged data sets.

  \code{head.join} and \code{tail.join} return a list, each component of
  which is the  result of a call to \code{\link{head}} or
  \code{\link{tail}} on each data set compont of the joined object.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{merge}}}
\examples{
## load the example data
data(swapdiat, swappH, rlgh)

## merge training and test set on columns
dat <- join(swapdiat, rlgh, verbose = TRUE)

## extract the merged data sets and convert to proportions
swapdiat <- dat[[1]] / 100
rlgh <- dat[[2]] / 100

## merge training and test set using left join
head(join(swapdiat, rlgh, verbose = TRUE, type = "left"))

## load the example data
data(ImbrieKipp, SumSST, V12.122)

## merge training and test set on columns
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- dat[[1]] / 100
V12.122 <- dat[[2]] / 100

## show just the first few lines of each data set
head(dat, n = 4)

## show just the last few lines of each data set
tail(dat, n = 4)

## merge training and test set using inner join
head(join(ImbrieKipp, V12.122, verbose = TRUE, type = "inner"))

## merge training and test set using outer join and replace
## NA with -99.9
head(join(ImbrieKipp, V12.122, verbose = TRUE, value = -99.9))
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
