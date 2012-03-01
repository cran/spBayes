\name{FBC07.dat}
\docType{data}

\alias{FBC07.dat}


\title{Synthetic multivariate data with spatial and non-spatial variance
structures}

\description{
The synthetic dataset describes a stationary and isotropic bivariate process. Please refer to the vignette Section 4.2 for specifics.

}

\usage{
data(FBC07.dat)
}
 
\format{
 A data frame of 250 rows and 4 columns. Columns 1 and 2 are
 coordinates and columns 3 and 4 are response variables.
}

\source{
  Finley A.O., S. Banerjee, and B.P. Carlin (2007) spBayes: R package for Univariate and Multivariate
  Hierarchical Point-referenced Spatial Models. Journal of Statistical Software.
}

\examples{

\dontrun{

data(FBC07.dat)
library(geoR)

max <- 40
bins <- 20
pts <- 1:150

vario.1 <- variog(coords=FBC07.dat[pts,1:2], data=FBC07.dat[pts,3],
          uvec=(seq(0, max, length=bins)))

vario.2 <- variog(coords=FBC07.dat[pts,1:2], data=FBC07.dat[pts,4],
                  uvec=(seq(0,max, length=bins)))

vario.fit.1 <-variofit(vario.1, ini.cov.pars=c(5.0, 1.0),
                       cov.model="exponential",
                       minimisation.function="nls",
                       weights="equal")

vario.fit.2 <-variofit(vario.2, ini.cov.pars=c(5.0, 10.0),
                       cov.model="exponential",
                       minimisation.function="nls",
                       weights="equal")

par(mfrow=c(1,2))

plot(vario.1$u, vario.1$v, axes=FALSE, type = "n",
     ylim=c(0,15), xlab="Distance", ylab="Semivariance")
points(vario.1$u, vario.1$v, pch=19, cex=0.5)
axis(1, seq(0,max,10))
axis(2, seq(0,15,5))
abline(h=vario.fit.1$nugget)
abline(h=vario.fit.1$cov.pars[1]+vario.fit.1$nugget)
abline(v=3/(1/vario.fit.1$cov.pars[2]))
lines(vario.fit.1)

plot(vario.2$u, vario.2$v, axes=FALSE, type = "n",
     ylim=c(0,15), xlab="Distance", ylab="")
points(vario.2$u, vario.2$v, pch=19, cex=0.5)
axis(1, seq(0,max,10))
abline(h=vario.fit.2$nugget)
abline(h=vario.fit.2$cov.pars[1]+vario.fit.2$nugget)
abline(v=3/(1/vario.fit.2$cov.pars[2]))
lines(vario.fit.2)
}
}

\keyword{datasets}

