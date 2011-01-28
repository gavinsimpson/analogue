\name{bootstrap.wa}
\alias{bootstrap.wa}
\alias{print.bootstrap.wa}
%\alias{print.bootstrap.wa}
\title{Bootstrap estimation and errors for WA models}
\description{
  Function to calculate bootstrap statistics for transfer function
  models such as bootstrap estimates, model RMSEP, sample specific
  errors for predictions and summary statistics such as bias and
  \eqn{R^2} between oberved and estimated environment.
}
\usage{
\method{bootstrap}{wa}(object, n.boot = 1000, verbose = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object of class \code{"wa"} for which bootstrap
    statistics are to be generated.}
  \item{n.boot}{numeric; the number of bootstrap samples to draw.}
  \item{verbose}{logical; should bootstrap progress be printed to the
    console?}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  See \code{\link{bootstrap.mat}} for further details. This method is
  not as feature packed as \code{\link{bootstrap.mat}} but can be used
  to evaluate the model performance of WA transfer function models.
}
\value{
  An object with the same components as \code{\link{predict.wa}}.
}
\references{

  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C. and ter Braak,
  C.J.F. (1990). Diatoms and pH reconstruction. \emph{Philosophical
    Transactions of the Royal Society of London; Series B}, \strong{327};
  263--278.

}
\author{Gavin L. Simpson }
\seealso{\code{wa}, \code{plot.wa}.}
\examples{
## continue the RLGH and SWAP example from ?wa
##example(wa)

## bootstrap the WA model
##boot.mod <- bootstrap(mod, n.boot = 100)

## performance statistics
##performance(boot.mod)

## Imbrie and Kipp
data(ImbrieKipp)
data(SumSST)
ik.wa <- wa(SumSST ~ ., data = ImbrieKipp, tol.dw = TRUE,
            min.tol = 2, small.tol = "min")
ik.wa

## compare actual tolerances to working values
with(ik.wa, rbind(tolerances, model.tol))

## bootstrap the WA model
ik.boot <- bootstrap(ik.wa, n.boot = 100)

## performance statistics
performance(ik.boot)
}
\keyword{multivariate}
\keyword{methods}
