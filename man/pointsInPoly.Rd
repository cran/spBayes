\name{pointsInPoly}
\alias{pointsInPoly}          
\title{Finds points in a polygon}

\description{
  Given a polygon and a set of points this function returns the subset
  of points that are within the polygon.
}

\usage{
 pointsInPoly(poly, points, ...)
}

\arguments{
  \item{poly}{an \eqn{n\times 2}{nx2} matrix of polygon vertices. Matrix
  columns correspond to vertices' x and y coordinates, respectively. }
  \item{points}{an \eqn{m\times 2}{mx2} matrix of points. Matrix
  columns correspond to points' x and y coordinates, respectively.}
  \item{...}{currently no additional arguments.}
}

\value{
  If points are found with the polygon, then a vector is returned with
  elements corresponding to the row indices of \code{points}, otherwise
  \code{NA} is returned.
}

\details{
It is assumed that the polygon is to be closed by joining the last vertex to the first vertex.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
}

\examples{

\dontrun{
##Example 1
points <- cbind(runif(1000, 0, 10),runif(1000, 0, 10))

poly <- cbind(c(1:9,8:1), c(1,2*(5:3),2,-1,17,9,8,2:9))

point.indx <- pointsInPoly(poly, points)

plot(points, pch=19, cex=0.5, xlab="x", ylab="y", col="red")
points(points[point.indx,], pch=19, cex=0.5, col="blue")
polygon(poly)

##Example 2
##a function to partition the domain
tiles <- function(points, x.cnt, y.cnt, tol = 1.0e-10){

  x.min <- min(points[,1])-tol
  x.max <- max(points[,1])+tol
  y.min <- min(points[,2])-tol
  y.max <- max(points[,2])+tol
  
  x.cnt <- x.cnt+1
  y.cnt <- y.cnt+1
  
  x <- seq(x.min, x.max, length.out=x.cnt)
  y <- seq(y.min, y.max, length.out=y.cnt)

  tile.list <- vector("list", (length(y)-1)*(length(x)-1))

  l <- 1
  for(i in 1:(length(y)-1)){
    for(j in 1:(length(x)-1)){
      tile.list[[l]] <- rbind(c(x[j], y[i]),
                              c(x[j+1], y[i]),
                              c(x[j+1], y[i+1]),
                              c(x[j], y[i+1]))
       l <- l+1
    }

  }

  tile.list
}

n <- 1000
points <- cbind(runif(n, 0, 10), runif(n, 0, 10))

grd <- tiles(points, x.cnt=10, y.cnt=10)

plot(points, pch=19, cex=0.5, xlab="x", ylab="y")

sum.points <- 0
for(i in 1:length(grd)){
  polygon(grd[[i]], border="red")
  
  point.indx <- pointsInPoly(grd[[i]], points)
  
  if(!is.na(point.indx[1])){
    sum.points <- length(point.indx)+sum.points
    
    text(mean(grd[[i]][,1]), mean(grd[[i]][,2]), length(point.indx), col="red")
  }
}
sum.points

}

}

\keyword{misc}
