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
\details{
Further information is available in the following vignettes:
\tabular{ll}{
\code{analogue_methods} \tab Analogue Methods in Palaeoecology (source, pdf)\cr
}
}
\section{Analogue methods}{
  \pkg{analogue} provides functions for analogue matching and the modern
  analogue technique (MAT) via \code{\link{analog}} and
  \code{\link{mat}}. A wide range of dissimilarity coefficients are
  available via the \code{\link{distance}} function.
}
\section{Transfer function methods}{
  Several approaches to fitting transfer function models are provided by
  \pkg{analogue}:
  \describe{
    \item{ \code{\link{wa}}: }{Simple and tolerance-downweighted weighted
      averaging with classical, inverse, and monotonic deshrinking.}
    \item{ \code{\link{mat}}: }{The modern analogue technique (MAT).}
    \item{ \code{\link{pcr}}: }{Principal components regression with
      ecologically meaningful transformations}
  }
}
\author{
Gavin L. Simpson, Jari Oksanen

Maintainer: Gavin L. Simpson <ucfagls@gmail.com>
}
\keyword{ package }
