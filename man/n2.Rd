\name{n2}
\alias{n2}
\alias{n2.default}

\title{Calculate Hill's N2 diversity measure}

\description{
  Hills N2 is a measure of species diversity, commonly referred to as
  "effective" diversity. If computed on the rows (samples) then the
  "effective" number of species in each row is returned, whereas, if
  computed on the columns (species) then the "effective" number of
  occurences of each species in the data set is returned.
}

\usage{
n2(x, \dots)

\method{n2}{default}(x, which = c("species", "sites"), \dots)
}


\arguments{
  \item{x}{matrix or data frame of species data}
  \item{which}{character; compute N2 on the rows (\code{"sites"}) or the
    columns (\code{"species"}) of \code{x}}
  \item{\dots}{arguments passed to other methods}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  A numeric vector of N2 values.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Gavin L. Simpson
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
data(swapdiat)
sppN2 <- n2(swapdiat, "species")
head(sppN2)
sampN2 <- n2(swapdiat, "sites")
head(sampN2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }
\keyword{ utilities }
