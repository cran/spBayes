\name{NYOzone.dat}
\docType{data}

\alias{NYOzone.dat}
\title{Observations of ozone concentration levels.}

\description{
  These data and subsequent description are drawn from the
  \pkg{spTimer} package (version 0.7). This data set contains values of daily 
  8-hour maximum average ozone concentrations (ppb; O3.8HRMAX), maximum  temperature (degree Celsius; cMAXTMP), wind speed
  (knots; WDSP), and 
  relative humidity (RM), obtained from 28 monitoring sites in New
  York, USA, between July 1 and August 31
  in 2006. Each row represents a station and columns hold
  consecutive daily values.
}

\usage{
data(NYOzone.dat)
}

\references{
  \pkg{spTimer} Bakar, K.S. and S.K. Sahu. \url{http://www.southampton.ac.uk/~sks/research/papers/spTimeRpaper.pdf}
  
  Sahu, S.K. and K.S. Bakar. (2012) A Comparison of Bayesian Models for Daily Ozone Concentration Levels. \emph{Statistical Methodology}, 9, 144--157.
}


\format{Columns for NYdata:
  \itemize{
    \item   1st col = Longitude
    \item   2nd col = Latitude
    \item   3rd col = X coordinates in UTM projection
    \item   4th col = Y coordinates in UTM projection
    \item   5th col = Ozone July 1 (O3.8HRMAX.1) 
    \item   6th col = Ozone July 2 (O3.8HRMAX.2)
    \item   ...
    \item   66th col = Ozone August 31 (O3.8HRMAX.62)
    \item   remaining columns organize cMAXTMP, WDSP, and RH
    identical to the 62 O3.8HRMAX measurements
  }
}
\keyword{datasets}
