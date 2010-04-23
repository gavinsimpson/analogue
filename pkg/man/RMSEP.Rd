\name{RMSEP}
\alias{RMSEP}
\alias{RMSEP.default}
\alias{RMSEP.mat}
\alias{RMSEP.bootstrap.mat}
\alias{RMSEP.bootstrap.wa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Root mean square error of prediction}
\description{
  Calculates or extracts the RMSEP from transfer function models.
}
\usage{
RMSEP(object, \dots)

\method{RMSEP}{mat}(object, k, weighted = FALSE,
          \dots)

\method{RMSEP}{bootstrap.mat}(object, type = c("birks1990", "standard"),
          \dots)

\method{RMSEP}{bootstrap.wa}(object, type = c("birks1990", "standard"),
          \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An R object.}
  \item{k}{numeric; the number of analogues to use in calculating the
    RMSEP. May be missing. If missing, \code{k} is extracted from the
    model using \code{\link{getK}}.}
  \item{weighted}{logical; Return the RMSEP for the weighted or
    unweighted model? The default is for an unweighted model.}
  \item{type}{The type of RMSEP to return/calculate. See Details,
          below.}
  \item{\dots}{Arguments passed to other methods.}
}
\details{
  There are two forms of RMSEP in common usage. Within palaeoecology,
  the RMSEP of Birks et al. (1990) is most familiar:
  \deqn{\mathrm{RMSEP} = \sqrt{s_1^2 + s_2^2}}{%
    RMSEP = sqrt(s[1]^2 + s[2]^2)}
  where where \eqn{s_1}{s[1]} is the standard deviation of the
  out-of-bag (OOB) residuals and \eqn{s_2}{s[2]} is the mean bias or the
  mean of the OOB residuals.

  In the wider statistical literature, the following form of RMSEP is
  more commonly used:
  \deqn{\mathrm{RMSEP} = \sqrt{\frac{\sum_{i=1}^n (y_i - \hat{y}_i)^2}{n}}}{%
    RMSEP = sqrt(mean((y[i] - yhat[i])^2))}
  where \eqn{y_i}{y[i]} are the observed values and \eqn{\hat{y}_i} the
  transfer function predictions/fitted values.

  The first form of RMSEP is returned by default or if \code{type =
    "birks1990"} is supplied. The latter form is returned if \code{type
    = "standard"} is supplied.

  The RMSEP for objects of class \code{"mat"} is a leave-one-out
  cross-validated RMSEP, and is calculated as for \code{type =
    "standard"}.
}
\value{
  A numeric vector of length 1 that is the RMSEP of \code{object}.
}
\references{
  
  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C. and ter Braak,
  C.J.F. (1990). Diatoms and pH reconstruction. \emph{Philosophical
    Transactions of the Royal Society of London; Series B}, \strong{327};
  263--278.
}
\author{Gavin L. Simpson}
\seealso{\code{\link{mat}}, \code{\link{bootstrap}}, \code{\link{wa}},
  \code{\link{bootstrap.wa}}.}
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

## fit the MAT model using the squared chord distance measure
(ik.mat <- mat(ImbrieKipp, SumSST, method = "chord"))

## Leave-one-out RMSEP for the MAT model
RMSEP(ik.mat)

## bootstrap training set
(ik.boot <- bootstrap(ik.mat, n.boot = 100))

## extract the Birks et al (1990) RMSEP
RMSEP(ik.boot)

## Calculate the alternative formulation
RMSEP(ik.boot, type = "standard")
}
\keyword{methods}
\keyword{utilities}
