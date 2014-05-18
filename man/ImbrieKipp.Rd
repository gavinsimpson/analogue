\name{ImbrieKipp}
\alias{ImbrieKipp}
\alias{SumSST}
\alias{WinSST}
\alias{Salinity}
\alias{V12.122}
\docType{data}
\title{Imbrie and Kipp foraminifera training set}
\description{
  The classic Imbrie and Kipp (1971) training set of counts on 27
  species of foraminifera from 61 ocean sediment core surface samples
  and associated measures of summer and winter sea-surface temperatures
  and salinity at each location.

  110 sediment cores samples from core V12-122 are also supplied in
  \code{V12.122}.
}
\usage{
data(SumSST)
data(WinSST)
data(Salinity)
data(V12.122)
}
\format{
  \code{ImbrieKipp} is a data frame with 61 observations on the
  following 27 species:
  \describe{
    \item{\code{O.univ}}{a numeric vector}
    \item{\code{G.cglob}}{a numeric vector}
    \item{\code{G.ruber}}{a numeric vector}
    \item{\code{G.tenel}}{a numeric vector}
    \item{\code{G.saccu}}{a numeric vector}
    \item{\code{G.rubes}}{a numeric vector}
    \item{\code{G.pacL}}{a numeric vector}
    \item{\code{G.pacR}}{a numeric vector}
    \item{\code{G.bullo}}{a numeric vector}
    \item{\code{G.falco}}{a numeric vector}
    \item{\code{G.calid}}{a numeric vector}
    \item{\code{G.aequi}}{a numeric vector}
    \item{\code{G.gluti}}{a numeric vector}
    \item{\code{G.duter}}{a numeric vector}
    \item{\code{G.infla}}{a numeric vector}
    \item{\code{G.trnL}}{a numeric vector}
    \item{\code{G.trnR}}{a numeric vector}
    \item{\code{G.crasf}}{a numeric vector}
    \item{\code{G.scitu}}{a numeric vector}
    \item{\code{G.mentu}}{a numeric vector}
    \item{\code{P.obliq}}{a numeric vector}
    \item{\code{C.nitid}}{a numeric vector}
    \item{\code{S.dehis}}{a numeric vector}
    \item{\code{G.digit}}{a numeric vector}
    \item{\code{Other}}{a numeric vector}
    \item{\code{G.quin}}{a numeric vector}
    \item{\code{G.hirsu}}{a numeric vector}
  }
  Summer and Winter sea-surface temperatures, and salinity values for
  the 61 sites in the Imbrie and Kipp training set (\code{ImbrieKipp}):
  \describe{
    \item{\code{SumSST}}{a numeric vector of summer sea-surface water
      temperatures. Values are in degrees C.}
    \item{\code{WinSST}}{a numeric vector of winter sea-surface water
      temperatures. Values are in degrees C.}
    \item{\code{Salinity}}{a numeric vector of sea water salinity
      values.}
  }
  \code{V12.122} is a data frame with 110 observations from core
  V12-122 on the following 28 species:
  \describe{
    \item{\code{O.univ}}{a numeric vector}
    \item{\code{G.cglob}}{a numeric vector}
    \item{\code{G.ruber}}{a numeric vector}
    \item{\code{G.tenel}}{a numeric vector}
    \item{\code{G.saccu}}{a numeric vector}
    \item{\code{G.rubes}}{a numeric vector}
    \item{\code{G.pacL}}{a numeric vector}
    \item{\code{G.pacR}}{a numeric vector}
    \item{\code{G.bullo}}{a numeric vector}
    \item{\code{G.falco}}{a numeric vector}
    \item{\code{G.calid}}{a numeric vector}
    \item{\code{G.aequi}}{a numeric vector}
    \item{\code{G.gluti}}{a numeric vector}
    \item{\code{G.duter}}{a numeric vector}
    \item{\code{G.infla}}{a numeric vector}
    \item{\code{G.trnL}}{a numeric vector}
    \item{\code{G.trnR}}{a numeric vector}
    \item{\code{G.crasf}}{a numeric vector}
    \item{\code{G.scitu}}{a numeric vector}
    \item{\code{G.mentu}}{a numeric vector}
    \item{\code{P.obliq}}{a numeric vector}
    \item{\code{C.nitid}}{a numeric vector}
    \item{\code{S.dehis}}{a numeric vector}
    \item{\code{G.digit}}{a numeric vector}
    \item{\code{G.hexag}}{a numeric vector}
    \item{\code{G.cglom}}{a numeric vector}
    \item{\code{cfH.pel}}{a numeric vector}
    \item{\code{Other}}{a numeric vector}
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Imbrie and Kipp (1971) TODO: Get the full reference.

  These data were provided in electronic format by Prof. H. John
  B. Birks.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(ImbrieKipp)
head(ImbrieKipp)

data(SumSST)
data(WinSST)
data(Salinity)

plot(cbind(SumSST, WinSST, Salinity))

data(V12.122)
head(V12.122)
}
\keyword{datasets}
