\name{tortula}
\alias{tortula}
\docType{data}
\title{Morphological data for ten taxa of the genus Tortula}
\description{
  These data are observations on a series of seven morphological
  variables for individuals in of the \emph{Tortula}
  sect. \emph{Rurales} De Not. (\emph{Pottiaceae, Musci}. 
}
\usage{
data(tortula)
}
\format{
  \code{tortula} is a data frame of seven morphological measurements on
  14 individuals from the genus \emph{Tortula}.
  \describe{
    \item{\code{Taxon}}{factor; the species of \emph{Tortula}}
    \item{\code{Hydroid}}{logical; presence of hydroid cells}
    \item{\code{LeafOutline}}{ordered; shape of the leaf outline}
    \item{\code{Denticulation}}{ordered; degree of denticulation}
    \item{\code{ApexShape}}{ordered; shape of the leaf apex}
    \item{\code{Length}}{numeric, leaf length}
    \item{\code{Diameter}}{numeric; leaf diameter}
    \item{\code{Papillae}}{numeric; number of papillae per cell}
  }
  The last three variables are the average of ten replicate samples from
  the same herbarium capsule.
}
\source{
  The data were presented in Podani (1999).
}
\references{
  Podani, J. (1999) Extending Gower's coefficient of similarity to
  ordinal characters. \emph{Taxon} \strong{48}, 331-340.
}
\examples{
data(tortula)
head(tortula)
str(tortula)
}
\keyword{datasets}
