\name{analogue-package}
\alias{analogue-package}
\alias{analogue}
\docType{package}
\title{
  Analogue and weighted averaging methods for palaeoecology
}
\description{
  \pkg{analogue} is a package for quantitative palaeoecology with a
  focus on analogue methods, transfer functions, and data handling and
  display.
}
\section{Analogue methods}{
  \pkg{analogue} provides functions for analogue matching and the modern
  analogue technique (MAT) via \code{\link{analog}} and
  \code{\link{mat}}. A wide range of dissimilarity coefficients are
  available via the \code{\link{distance}} function.

  Additional analysis of modern and no-analogue problems is facilitated
  via a range of functions implementing many methods from the
  literature. In particular, the receiver operating characteric (ROC)
  curves method of Gavin et al (2003) is available in \code{\link{roc}}
  and a related method employing direct logistic regression modelling
  (Simpson & Birks, 2012) is available in \code{\link{logitreg}}.
}
\section{Transfer function methods}{
  Several approaches to fitting transfer function models are provided by
  \pkg{analogue}:
  \describe{
    \item{ \code{\link{wa}}: }{Simple and tolerance-downweighted weighted
      averaging with classical, inverse, and monotonic spline
      deshrinking.}
    \item{ \code{\link{mat}}: }{The modern analogue technique (MAT).}
    \item{ \code{\link{pcr}}: }{Principal components regression with
      ecologically meaningful transformations}
  }

  A range of functions for working with and exploring training sets and
  palaeoenvironmental reconstructions is also included in
  \pkg{analogue}. These include

  \describe{
    \item{ \code{\link{crossval}} }{leave-one-out, repeated k-fold, and
    bootstrap cross-validation methods.}
    \item{ \code{\link{compare}} }{compare properties of taxa or other
      proxies across modern and fossil data sets.}
    \item{ \code{\link{evenSample}} }{are training set samples evenly
      distributed along the gradient of interest?}
    \item{ \code{\link{splitSample}} }{splits a gradient into a set of
      bins or chunks and samples evenly from within each chunk to create
      a representative test set for cross-validation.}
    \item{ \code{\link{timetrack}} }{overlays a fossil or secondary data
      set on to an (constrained) ordination of a modern or reference
      data set.}
    \item{ \code{\link{weightedCor}} }{implements the weighted
      correlation test of a Weighted Averaging reconstruction as proposed
      by Telford & Birks (2011).}
  }
}
\section{Utilities}{
  \pkg{analogue} provides a range of utilities for working with palaeo
  data.

  \describe{
    \item{ \code{\link{tran}} }{a range of transformations applicable to
      or commonly used with palaeo data.}
    \item{ \code{\link{Stratiplot}} }{draws stratigraphic diagrams using
      the Lattice package.}
    \item{ \code{\link{join}} }{merging of modern/training set and
      fossil data sets.}
    \item{ \code{\link{chooseTaxa}} }{selects taxa that meet certain
      abundance and occurrence criteria.}
  }
}
\section{Documentation}{
  A full tutorial and worked example for the main features of analogue
  matching and MAT is avilable in the vignette

  \describe{
    \item{\code{analogue_methods}}{Analogue Methods in Palaeoecology}
  }
}
\author{
Gavin L. Simpson, Jari Oksanen

Maintainer: Gavin L. Simpson <ucfagls@gmail.com>
}
\references{
  
  Gavin, D.G., Oswald, W.W., Wahl, E.R. and Williams, J.W. (2003) A
  statistical approach to evaluating distance metrics and analog
  assignments for pollen records. \emph{Quaternary Research}
  \strong{60}, 356--367.

  Simpson, G.L. & Birks H.J.B. (2012) Statistical Learning in
  Palaeolimnology. In Birks, H.J.B, Lotter, A.F. Juggins S., and Smol,
  J.P. (Eds) \emph{Tracking Environmental Change Using Lake Sediments,
    Volume 5: Data Handling and Numerical Techniques}. Springer,
  Dordrecht.
  
  Telford R.J. and Birks, H.J.B. (2011) A novel method for assessing the
  statistical significance of quantitative reconstructions inferred from
  biotic assemblages. \emph{Quanternary Science Reviews}
  \strong{30}:1272-1278.
}
\keyword{ package }
