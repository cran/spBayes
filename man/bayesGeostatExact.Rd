\name{bayesGeostatExact}
\alias{bayesGeostatExact}          
\title{Simple Bayesian spatial linear model with fixed semivariogram parameters }

\description{
  Given a observation coordinates and fixed semivariogram
  parameters the \code{bayesGeostatExact} function fits a
  simple Bayesian spatial linear model.  
}

\usage{
  bayesGeostatExact(formula, data = parent.frame(), n.samples,
                     beta.prior.mean, beta.prior.precision,
                     coords, cov.model="exponential", phi, nu, alpha,
                     sigma.sq.prior.shape, sigma.sq.prior.rate,
                     sp.effects=TRUE, verbose=TRUE, ...)
}

\arguments{
  \item{formula}{for a univariate model, this is a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spLM} is called.  }
  \item{n.samples}{the number of posterior samples to collect. }
  \item{beta.prior.mean}{\eqn{\beta}{beta} multivariate normal mean vector hyperprior. }
  \item{beta.prior.precision}{\eqn{\beta}{beta} multivariate normal precision
    matrix hyperprior. }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{cov.model}{a quoted key word that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{phi}{the fixed value of the spatial decay. }
  \item{nu}{if \code{cov.model} is \code{"matern"} then the fixed value
    of the spatial process smoothness must be specified. }
  \item{alpha}{the fixed value of the ratio between the nugget
    \eqn{\tau^2}{tau.sq} and partial-sill \eqn{\sigma^2}{sigma.sq}
    parameters from the specified \code{cov.model}. }  
  \item{sigma.sq.prior.shape}{\eqn{\sigma^2}{sigma.sq} (i.e., partial-sill) inverse-Gamma shape
    hyperprior. }
  \item{sigma.sq.prior.rate}{\eqn{\sigma^2}{sigma.sq} (i.e., partial-sill) inverse-Gamma 1/scale
    hyperprior. }
  \item{sp.effects}{a logical value indicating if spatial random effects
    should be recovered. }
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{...}{currently no additional arguments.  }
}
 
\value{
    An object of class \code{bayesGeostatExact}, which is a list with the following tags:
  \item{p.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{sp.effects}{a matrix that holds samples from the posterior
    distribution of the spatial random effects. The rows of this matrix
    correspond to the \eqn{n}{n} point observations and the columns are the
    posterior samples. }
  \item{args}{a list with the initial function arguments.  }
}

%\references{
  %\url{http://blue.fr.umn.edu/spatialBayes}
%}
  
\author{
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\examples{
\dontrun{

data(FBC07.dat)
Y <- FBC07.dat[1:150,"Y.2"]
coords <- as.matrix(FBC07.dat[1:150,c("coord.X", "coord.Y")])

n.samples <- 500
n = length(Y)
p = 1

phi <- 0.15
nu <- 0.5

beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

alpha <- 5/5

sigma.sq.prior.shape <- 2.0
sigma.sq.prior.rate <- 5.0

##############################
##Simple linear model with
##the default exponential
##spatial decay function
##############################
set.seed(1)
m.1 <- bayesGeostatExact(Y~1, n.samples=n.samples,
                          beta.prior.mean=beta.prior.mean,
                          beta.prior.precision=beta.prior.precision,
                          coords=coords, phi=phi, alpha=alpha,
                          sigma.sq.prior.shape=sigma.sq.prior.shape,
                          sigma.sq.prior.rate=sigma.sq.prior.rate)



print(summary(m.1$p.samples))

##Requires MBA package to
##make surfaces
library(MBA)
par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, Y), no.X=100, no.Y=100, extend=T)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

w.hat <- rowMeans(m.1$sp.effects)
w.surf <-
  mba.surf(cbind(coords, w.hat), no.X=100, no.Y=100, extend=T)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Estimated random effects")
points(coords)
contour(w.surf, add=T)


##############################
##Simple linear model with
##the matern spatial decay
##function. Note, nu=0.5 so
##should produce the same
##estimates as m.1
##############################
set.seed(1)
m.2 <- bayesGeostatExact(Y~1, n.samples=n.samples,
                          beta.prior.mean=beta.prior.mean,
                          beta.prior.precision=beta.prior.precision,
                          coords=coords, cov.model="matern",
                          phi=phi, nu=nu, alpha=alpha,
                          sigma.sq.prior.shape=sigma.sq.prior.shape,
                          sigma.sq.prior.rate=sigma.sq.prior.rate)

print(summary(m.2$p.samples))

##############################
##This time with the
##spherical just for fun
##############################
m.3 <- bayesGeostatExact(Y~1, n.samples=n.samples,
                          beta.prior.mean=beta.prior.mean,
                          beta.prior.precision=beta.prior.precision,
                          coords=coords, cov.model="spherical",
                          phi=phi, alpha=alpha,
                          sigma.sq.prior.shape=sigma.sq.prior.shape,
                          sigma.sq.prior.rate=sigma.sq.prior.rate)

print(summary(m.3$p.samples))

##############################
##Another example but this
##time with covariates
##############################
data(FORMGMT.dat)

n = nrow(FORMGMT.dat)
p = 5 ##an intercept an four covariates

n.samples <- 50

phi <- 0.0012

coords <- cbind(FORMGMT.dat$Longi, FORMGMT.dat$Lat)
coords <- coords*(pi/180)*6378

beta.prior.mean <- rep(0, times=p)
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

alpha <- 1/1.5

sigma.sq.prior.shape <- 2.0
sigma.sq.prior.rate <- 10.0

m.4 <-
  bayesGeostatExact(Y~X1+X2+X3+X4, data=FORMGMT.dat, n.samples=n.samples,
                     beta.prior.mean=beta.prior.mean,
                     beta.prior.precision=beta.prior.precision,
                     coords=coords, phi=phi, alpha=alpha,
                     sigma.sq.prior.shape=sigma.sq.prior.shape,
                     sigma.sq.prior.rate=sigma.sq.prior.rate)

print(summary(m.4$p.samples))



##Requires MBA package to
##make surfaces
library(MBA)
par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, resid(lm(Y~X1+X2+X3+X4, data=FORMGMT.dat))),
                 no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

w.hat <- rowMeans(m.4$sp.effects)
w.surf <-
  mba.surf(cbind(coords, w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Estimated random effects")
contour(w.surf, add=T)
points(coords, pch=1, cex=1)


}
}

\keyword{misc}
