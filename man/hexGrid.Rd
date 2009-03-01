\name{hexGrid}
\alias{hexGrid}
\title{Generates a hexagon tessellation within a bounding box }
\description{
  This function generates a hexagon tessellation within a bounding box.
}
\usage{
 hexGrid(nodes.per.layer, b.box, ...)
}
\arguments{
  \item{nodes.per.layer}{ the number of hexagons from \code{min.x} to \code{max.x} at \code{min.y}. }
  \item{b.box}{ the domain's bounding box defined with a vector of length
    four and elements ordered \code{min.x}, \code{min.y}, \code{max.x}, \code{max.y}. }
  \item{\dots}{ currently no additional arguments. }
}

\value{
  A list with the following tags:
  \item{hx.hy}{ a vector of length two that holds the height and width,
    respectively, between hexagon centroids. }
  \item{layers}{ the number of hexagon layers within the domain. }
  \item{nodes.per.layer}{ the number of hexagons from \code{min.x} to \code{max.x} at \code{min.y}. }
  \item{hex.centroids}{ an \eqn{n\times 2}{nx2} matrix with rows
    corresponding the the \eqn{x}{x} and \eqn{y}{y} coordinates of the
    \eqn{n}{n} hexagons within the domain. }
  \item{hex.polygons}{ a list with \eqn{6\times 2}{6x2} matrices which
    hold the hexagons' \eqn{x}{x} and \eqn{y}{y} vertices coordinates.}
}

\author{Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}}

\examples{
\dontrun{

##Define the bounding box and make the call
b.box <- c(0,0,10,10)

out <- hexGrid(20, b.box)

##Plot using lapply
par(mfrow=c(1,2))
    
plot(out$hex.centroids, pch=19, cex=0.5,
     ylab="Northing", xlab="Easting")
lapply(out$hex.polygons, polygon, col="blue")

##Now color hexagons based on value
my.col.ramp <- function(z){
  zlim <- range(z)
  zlen <- zlim[2]-zlim[1]+1
  colorlut <- heat.colors(as.integer(zlen))
  col <- colorlut[z-zlim[1]+1]
  col
}

n <- nrow(out$hex.centroids)
col <- my.col.ramp(rnorm(n))

plot(out$hex.centroids, typ="n", ylab="", xlab="Easting")
for(i in 1:n){
  polygon(out$hex.polygons[[i]], col=col[i], border=col[i])
}

}
  }


\keyword{misc}
