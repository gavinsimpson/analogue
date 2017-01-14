\name{Pollen}
\alias{Pollen}
\alias{Biome}
\alias{Climate}
\alias{Location}
\docType{data}
\title{North American Modern Pollen Database}
\description{
  A database of modern pollen samples from a network of sites from North
  America and Greenland, compiled by Whitmore et al. (2005). Associated
  climatic and vegetation data are also record. The version of the NAMPD
  included here is latest version 1-7.3 (February 2013), as of January
  2016.
}
\usage{
data(Pollen)

data(Climate)

data(Biome)

data(Location)
}
\format{
  For \code{Pollen}, a data frame of 4833 samples and 135 columns (the
  unique identifier and 134 pollen taxa).
  
  For \code{Biome}, a data frame of 4833 samples on, currently, a single
  vegetation variable (plus the unique identifier):
  \describe{
    \item{\code{ID2}}{Unique, sequential number assigned by NAMPD.}
    \item{\code{Fedorova}}{Factor; Vegetation type (Biome) See Whitmore
      et al. (2005), Figs 3 & 4. Reclassified biomes from Fedorova et al
      (1994).}
  }

  For \code{Location}, a data frame of the latitude and longitude
  locational data for 4833 samples.
  \describe{
    \item{\code{ID2}}{Unique, sequential number assigned by NAMPD.}
    \item{\code{Latitude}}{Latitude of the sampling location in decimal
      degrees.}
    \item{\code{Longitude}}{Longitude of the sampling location in decimal
      degrees.}
  }

  For \code{Climate}, a data frame with 4833 observations on the
  following 32 variables.
  \describe{
    \item{\code{ID2}}{Unique, sequential number assigned by NAMPD.}
    \item{\code{t[jan:dec]}}{numeric vectors; Mean monthly temperatures
      for the indicated month. Degrees C.}
    \item{\code{p[jan:dec]}}{numeric vectors; Mean total monthly
      precipitation (mm) for the indicated month.}
    \item{\code{tave}}{numeric; annual average temperature (Degrees C)}
    \item{\code{tmax}}{numeric; maximum temperature (in Degrees C)
      observed over the period of record.}
    \item{\code{tmin}}{numeric; minimum temperature (in Degrees C)
      observed over the period of record.}
    \item{\code{gdd0}}{numeric; Growing degree days computed using a
      base of 0 degrees C.}
    \item{\code{gdd5}}{numeric; Growing degree days computed using a
      base of 5 degrees C.}
    \item{\code{mtco}}{numeric; mean temperature of the coldest month.}
    \item{\code{mtwa}}{numeric; mean temperature of the warmest month.}
    \item{\code{annp}}{numeric; mean annual total precipitation (mm).}
  }
}
\details{
  These datasets were extracted from Version 1.7 of the North American
  Modern Pollen Database.

  All pollen species were included, however, only the Vegetation type
  (Biome) field of the AVHRR data and selected Climatic variables were
  extracted. Requests for additional variables to be included in the
  versions of the data included in the package should me sent to the
  package maintainer.
}
\section{Warning}{
  Note that the data for the pollen species are a mixture of types. The
  \code{DataForm} variable contains information on the type of data
  included for each site. The codes are:
  \describe{
    \item{RC}{raw counts}
    \item{RP}{raw counts expressed as percentages}
    \item{DC}{digitised counts}
    \item{DP}{digitised counts expressed as percentages}
    \item{PM}{counts in permille}
  }
  This value is not known for all samples.
}
\source{
  The database is archived electronically at two sites:

  \url{http://www.lpc.uottawa.ca/data/modern/index.html}

  \url{http://www.geography.wisc.edu/faculty/williams/web/data.htm}
}
\references{
  Whitmore, J., Gajewski, K., Sawada, M., Williams, J.W., Shuman, B.,
  Bartlein, P.J., Minckley, T., Viau, A.E., Webb III, T., Anderson,
  P.M., and Brubaker L.B., 2005. North American and Greenland modern
  pollen data for multi-scale paleoecological and paleoclimatic
  applications. \emph{Quaternary Science Reviews}
  \strong{24}:1828--1848.
  
  Williams, J.W., Shuman, B., Bartlein, P.J., Whitmore, J., Gajewski,
  K., Sawada, M., Minckley, T., Shafer, S., Viau, A.E., Webb, III, T.,
  Anderson, P.M., Brubaker, L.B., Whitlock, C. and Davis, O.K., 2006.
  \emph{An Atlas of Pollen-Vegetation-Climate Relationships for the
    United States and Canada}. American Association of Stratigraphic
  Palynologists Foundation, Dallas, TX, 293p.

  Williams, J.W. and Shuman, B., 2008. Obtaining accurate and precise
  environmental reconstructions from the modern analog technique and
  North American surface pollen dataset. \emph{Quaternary Science
    Reviews}, \strong{27}: 669--687.

  Fedorova, I.T., Volkova, Y.A., Varlyguin, E., 1994. \emph{World
    vegetation cover. Digital raster data on a 30-minute cartesian
    orthonormal geodetic (lat/long) 1080x2160 grid.} In: Global
  Ecosystems Database Version 2.0. USDOC/NOAA National Geophysical
  Data Center, Bould, CO.
}
\examples{
data(Pollen)

data(Climate)

data(Biome)

data(Location)
}
\keyword{datasets}
