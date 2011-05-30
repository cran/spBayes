\name{spLM}
\alias{spLM}
\title{Function for fitting univariate Bayesian spatial regression models}

\description{
  The function \code{spLM} fits Gaussian univariate Bayesian
  spatial regression models. Given a set of knots, \code{spLM} will also
  fit a \emph{predictive process} model (see references below).
}

\usage{
spLM(formula, data = parent.frame(), coords, knots,
      starting, sp.tuning, priors, cov.model,
      modified.pp = TRUE, n.samples, sub.samples,
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spLM} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing)
    or a vector of length two or three with the first and second elements recording the
    number of columns and rows in the desired knot grid. The third,
    optional, element sets the offset of the outermost knots from the
    extent of the \code{coords} extent. } 
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}.  The value portion of each
    of each tag is the parameter's starting value. }
  \item{sp.tuning}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}. The value portion of each
    of each tag defines the variance of the Metropolis normal proposal distribution.}
  \item{modified.pp}{a logical value indicating if the \emph{modified
    predictive process} should be used (see references below for
    details). Note, if a predictive process model is not used (i.e., \code{knots} is not specified) then
    this argument is ignored. }
  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{sigma.sq.ig},
    \code{tau.sq.ig}, \code{phi.unif}, and \code{nu.unif} (\code{Beta} priors are
    assumed \emph{flat}). Variance parameters, \code{simga.sq} and
    \code{tau.sq}, are assumed to follow an
    inverse-Gamma distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Gamma are
    passed as a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length two with the first
    and second elements corresponding to the lower and upper support,
    respectively. }
  \item{cov.model}{a quoted key word that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{n.samples}{the number of MCMC iterations. }
  \item{sub.samples}{a vector of length 3 that specifies \emph{start},
  \emph{end}, and \emph{thin}, respectively, of the MCMC samples. The
  default is \code{c(1, n.samples, 1)} (i.e., all samples). }
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report Metropolis acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{spLM}, which is a list with the following
  tags:
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{knot.coords}{the \eqn{m \times 2}{m x 2} matrix as specified by \code{knots}.  }
  \item{p.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis sampling
    acceptance rate. }
  \item{sp.effects}{a matrix that holds samples from the posterior
    distribution of the spatial random effects. The rows of this matrix
    correspond to the \eqn{n}{n} point observations and the columns are the
    posterior samples. }
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\seealso{
  \code{\link{spGGT}}, \code{\link{spMvLM}}
}

\references{
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.

  Finley, A.O., S. Banerjee, P. Waldmann, and T. Ericsson. (2008).
  Hierarchical spatial modeling of additive and dominance genetic
  variance for large spatial trial datasets. \emph{Biometrics}. DOI: 10.1111/j.1541-0420.2008.01115.x

  Finley, A.O,. H. Sang, S. Banerjee, and A.E. Gelfand. (2008). Improving the performance of predictive process modeling for
large datasets. \emph{Computational Statistics and Data Analysis}, DOI: 10.1016/j.csda.2008.09.008    
  
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
w <- rf.n200.dat$w

##############################
##Simple spatial regression
##############################
m.1 <- spLM(Y~1, coords=coords,
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100)

print(summary(m.1$p.samples))
plot(m.1$p.samples)

##Requires MBA package to
##make surfaces
library(MBA)
par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, Y), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

w.hat <- rowMeans(m.1$sp.effects)
w.surf <-
  mba.surf(cbind(coords, w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Estimated random effects")
points(coords)
points(m.1$knot.coords, pch=19, cex=1)
contour(w.surf, add=T)

##############################
##Predictive process
##############################
##Use some more observations
data(rf.n500.dat)

Y <- rf.n500.dat$Y
coords <- as.matrix(rf.n500.dat[,c("x.coords","y.coords")])
w <- rf.n500.dat$w

m.2 <- spLM(Y~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             modified.pp=FALSE, n.samples=2000, verbose=TRUE, n.report=100)

print(summary(m.2$p.samples))
plot(m.2$p.samples)

##Requires MBA package to
##make surfaces
library(MBA)
par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

w.hat <- rowMeans(m.2$sp.effects)
w.surf <-
  mba.surf(cbind(coords, w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Estimated random effects")
contour(w.surf, add=T)
points(coords, pch=1, cex=1)
points(m.2$knot.coords, pch=19, cex=1)
legend(1.5,2.5, legend=c("Obs.", "Knots"), pch=c(1,19), bg="white")

##############################
##Modified predictive process
##############################
m.3 <- spLM(Y~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=2000, verbose=TRUE, n.report=100)

print(summary(m.3$p.samples))
plot(m.3$p.samples)

##Requires MBA package to
##make surfaces
library(MBA)
par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

w.hat <- rowMeans(m.3$sp.effects)
w.surf <-
  mba.surf(cbind(coords, w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Estimated random effects")
contour(w.surf, add=T)
points(coords, pch=1, cex=1)
points(m.3$knot.coords, pch=19, cex=1)
legend(1.5,2.5, legend=c("Obs.", "Knots"), pch=c(1,19), bg="white")


}
}

\keyword{misc}
