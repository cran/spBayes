\name{Zurich.dat}
\docType{data}

\alias{Zurich.dat}


\title{Zurichberg Forest inventory data}

\description{
  
  Inventory data of the Zurichberg Forest,
  Switzerland (see Mandallaz 2008 for details). These data are provided
  with the kind authorization of the Forest Service of the Caton of
  Zurich.
  
  This dataset holds the coordinates for all trees in the Zurichberg
  Forest. Species (SPP), basal area (BAREA) diameter at breast height (DBH),
  and volume (VOL) are recorded for each tree. See species codes below. 

}

\usage{
data(Zurich.dat)
}
 
\format{
 A data frame containing 4954 rows and 6 columns.
}

\examples{
\dontrun{
data(Zurich.dat)

coords <- Zurich.dat[,c("X_TREE", "Y_TREE")]

spp.name <- c("beech","maple","ash","other broadleaves",
              "spruce","silver fir", "larch", "other coniferous")

spp.col <- c("yellow","red","orange","pink",
             "green","dark green","black","gray")
                 
plot(coords, col=spp.col[Zurich.dat$SPP+1],
     pch=19, cex=0.5, ylab="Northing", xlab="Easting")

legend.coords <- c(23,240)

legend(legend.coords, pch=19, legend=spp.name,
       col=spp.col, bty="n")


}
}

\keyword{datasets}
