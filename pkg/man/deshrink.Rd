\name{deshrink}
\Rdversion{1.1}
\alias{deshrink}
\alias{deshrinkPred}
\title{
Deshrinking techniques for WA transfer functions
}
\description{
In Weighted Averaging models averages are taken twice and thus WA
estimates shrink towards the training set mean and need to be
deshrunk.\code{deshrink} performs this deshrinking using several
techniques, whilst \code{deshrinkPred} will deshrink WA estimates for
new samples given a set of deshrinking coefficients.
}
\usage{
deshrink(env, wa.env, type = c("inverse", "classical",
                               "expanded", "none"))

deshrinkPred(x, coef, type = c("inverse", "classical",
                               "expanded", "none"))
}
\arguments{
  \item{env}{numeric; original environmental values.}
  \item{wa.env}{numeric; initial weighted average estimates.}
  \item{type}{character; the type of deshrinking. One of
    \code{"inverse"}, \code{"classical"}, \code{"expand"},
    \code{"none"}.}
  \item{x}{numeric; estimates to be deshrunk.}
  \item{coef}{numeric; deshrinking coefficients to use. Currently needs
    to be a vector of length 2. These should be supplied in the order
    \eqn{\beta_0,\beta_1}{beta[0], beta[1]}.}
}
%%\details{
%% FIX ME - give details on the approach used
%%}
\value{
  For \code{deshrinkPred} a numeric vector of deshrunk estimates.
  
  For an object of class \code{"deshrink"}, inheriting from class
  \code{"list"}, with two components. The type of deshrinking performed
  is stroed within attribute \code{"type"}. The componets of the
  returned object are:
  \item{coefficients}{The deshrinking coefficients used.}
  \item{env}{The deshrunk WA estimates.}
}
\references{
  Birks, H.J.B. (1995) Quantitative environmental reconstructions. In
  \emph{Statistical modelling of Quaternary science data} (eds.~D.Maddy
  \& J.S. Brew). Quaternary Research Association technical guide
  5. Quaternary Research Association, Cambridge.
}
\author{
  Gavin L. Simpson \& Jari Oksanen
}
%%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{\code{deshrinkPred}, does not currently check that
  the correct coefficients have been supplied in the correct order.}

\seealso{
\code{\link{wa}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
