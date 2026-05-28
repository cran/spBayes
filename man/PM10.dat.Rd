\name{PM10.dat}
\docType{data}

\alias{PM10.dat}
\title{Observed and modeled PM10 concentrations across Europe}

\description{

The \code{PM10.dat} data frame is a subset of data analyzed in Hamm
et al. (2015) and Datta et al. (2016). Data comprise April 6, 2010
square root transformed PM10 measurements across
central Europe with corresponding output from the LOTOS-EUROS Schaap et
al. (2008) chemistry transport model (CTM). CTM data may differ
slightly from that considered in the studies noted above due to LOTOS-EUROS CTM
code updates. A \code{NA} value is given
at CTM output locations were PM10 is not
observed. Point coordinates are in
"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs".
}

\usage{
data(PM10.dat)
}

\references{
Datta A., S. Banerjee, A.O. Finley, N. Hamm, and M. Schaap (2016). Nonseparable dynamic nearest neighbor Gaussian process models for large spatio-temporal data with an application to particulate matter analysis. \emph{Annals of Applied Statistics}, 10(3), 1286--1316. ISSN 1932-6157. doi:10.1214/16-AOAS931.
  
Hamm N. A.O. Finley, M. Schaap, A. Stein (2015). A Spatially Varying
Coefficient Model for Mapping PM10 Air Quality at the European scale. \emph{Atmospheric Environment}, 102, 393--405.

  Schaap M., R.M.A Timmermans, M. Roemer, G.A.C. Boersen, P. Builtjes, F. Sauter, G. Velders, J. Beck (2008). The LOTOS-EUROS Model: Description, Validation and Latest Developments. \emph{International Journal of Environment and Pollution}, 32(2), 270--290.
}


\format{Columns for PM10.dat:
  \itemize{
    \item x.coord = x coordinate (see projection information in the description)
    \item y.coord = y coordinate (see projection information in the description)
    \item pm10.obs = square root transformed PM10 measurements at
  monitoring stations (\code{NA} means there is not a station at the
  given location)
    \item pm10.ctm = square root transformed PM10 from CTM
  }
}
\keyword{datasets}
