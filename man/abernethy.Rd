\name{abernethy}
\alias{abernethy}
\docType{data}
\title{Abernethy Forest Pollen Sequence}
\description{
  The classic pollen data set from Abernethy Forest in the Scottish
  highlands, UK. The data originate from the work of Hilary Birks and
  Rolf Mathewes (1978) and have been analysed in several texts on
  quantitative numerical palaeoecology.

  The data set consists of 36 pollen taxa from 49 levels, with two
  additional variables; \code{Age}, the age of each sample, and
  \code{Depth} the depth (in cm) below the surface of the peat sequence
  from which the core was taken.
}
\usage{
data(abernethy)
}
\format{
  \code{abernethy} is a data frame with 49 samples on 36 species plus
  sample Age and Depth (in cm).
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  These data were provided in electronic format by Prof. H. John
  B. Birks. The original source is Birks and Mathewes (1978).
}
\references{
  Birks, H.H. and Mathewes, R.W. (1978) Studies in the vegetational
  history of Scotland. \emph{New Phytologist} \strong{80}, 455-484.
}
\examples{
data(abernethy)
head(abernethy)

(plt <- Stratiplot(Age ~ . - Depth,
                   data = chooseTaxa(abernethy, n.occ = 5, max.abun = 10),
                   type = "poly"))
}
\keyword{datasets}
