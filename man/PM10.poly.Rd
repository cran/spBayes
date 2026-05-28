\name{PM10.poly}
\docType{data}

\alias{PM10.poly}
\title{European countries used in PM10.dat}

\description{
European countries corresponding to \code{PM10.dat} locations and used in Hamm et al. (2015) and Datta et al. (2016). Polygon projection is "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs".
}

\usage{
data(PM10.poly)
}

\references{
Datta A., S. Banerjee, A.O. Finley, N. Hamm, and M. Schaap (2016). Nonseparable dynamic nearest neighbor Gaussian process models for large spatio-temporal data with an application to particulate matter analysis. \emph{Annals of Applied Statistics}, 10(3), 1286--1316. ISSN 1932-6157. doi:10.1214/16-AOAS931.
  
Hamm N. A.O. Finley, M. Schaap, A. Stein (2015). A Spatially Varying
Coefficient Model for Mapping PM10 Air Quality at the European scale. \emph{Atmospheric Environment}, 102, 393--405.

}

\format{
  List of polygons. See example below to convert to a \code{SpatialPolygons} object.


}

\examples{

\dontrun{

library(sp)

prj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs"

pm10.poly <- SpatialPolygons(PM10.poly, pO = 1:length(PM10.poly), proj4string=CRS(prj))
}

}

\keyword{datasets}
