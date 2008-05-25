\name{sp.lm}
\alias{sp.lm}
\title{Function for fitting univariate Bayesian spatial regression models}

\description{
  The function \code{sp.lm} fits Gaussian univariate stationary Bayesian
  spatial regression models. Given a set of knots, \code{sp.lm} fits a
  \emph{predictive process} model (see references below).
}

\usage{
sp.lm(formula, data = parent.frame(), coords, knots,
      fixed, starting, sp.tuning, priors, cov.model,
      sp.effects=TRUE, modified.pp = TRUE, n.samples, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{for a univariate model, this is a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{sp.lm} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{Either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing)
    or a vector of length two or three with the first and second elements recoding the
    number of columns and rows in the desired knot grid. The third,
    optional, element sets the offset of the outermost knots from the
    extent of the \code{coords} extent. }
  \item{fixed}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}.  The value portion of each
    of each tag is value at which the parameter will be fixed. } 
  \item{starting}{a list with each tag corresponding to a
    parameter name (not already present in the \code{fixed} list). Valid list tags are \code{beta}, \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}.  The value portion of each
    of each tag is the parameter's starting value. }
  \item{sp.tuning}{a list with each tag corresponding to a variance,
    spatial range, or smoothness parameter names (not already present in the \code{fixed} list). Valid list tags are \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}.  The value portion of each
    of each tag defines the step size of the proposal used in the
    Metropolis-Hastings sampling. }
  \item{modified.pp}{a logical value indicating if the \emph{modified
    predictive process} should be used (see references below for
    details). Note, if a predictive process model is not used (i.e., \code{knots} is not specified) then
    this argument is ignored. }
  \item{priors}{a list with each tag corresponding to a
    parameter name (not already present in the \code{fixed} list). Valid list tags are \code{sigma.sq.ig},
    \code{tau.sq.ig}, \code{phi.unif}, and \code{nu.unif} (\code{Beta} priors are
    assumed \emph{flat}). Variance parameters, \code{simga.sq} and
    \code{tau.sq}, are assumed to follow an
    inverse-Gamma distribution, whereas the spatial range \code{phi}
    and smoothness \code{nu} parameters is assumed to follow a Uniform distribution. The
    hyperparameters of the inverse-Gamma are
    passed as a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperpriors
    of the Uniform are also passed as a vector of length two with the first
    and second elements corresponding to the lower and upper support,
    respectively. }
  \item{cov.model}{a quoted key word that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{sp.effects}{a logical value indicating if spatial random effects should be recovered. }
  \item{n.samples}{the number of MCMC iterations. }
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report MH acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{sp.lm}, which is a list with the following
  tags:
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{knot.coords}{the \eqn{m \times 2}{m x 2} matrix as specified by \code{knots}.  }
  \item{p.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis-Hastings sampling
    acceptance rate. }
  \item{sp.effects}{a matrix that holds samples from the posterior
    distribution of the spatial random effects. The rows of this matrix
    correspond to the \eqn{n}{n} point observations and the columns are the
    posterior samples. }
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\seealso{
  \code{\link{ggt.sp}}
}

\references{
  Banerjee, S., Gelfand, A.E., Finley, A.O., and Sang, H. (In press). Gaussian predictive process models for large spatial datasets. \emph{Journal of the Royal Statistical Society Series B}.

  Finley, A.O., S. Banerjee, P. Waldmann, and T. Ericsson. (2007)
  Hierarchical spatial modeling of additive and dominance genetic
  variance for large spatial trial datasets. Research Report, Division
  of Biostatistics, University of Minnesota, 2007. \emph{Biometrics}, In
  press.

  Finley, A.O,. H. Sang, S. Banerjee, and A.E. Gelfand. (2008) Predictive process models for large multivariate spatial datasets. Research Report, Division
  of Biostatistics, University of Minnesota, 2008. \emph{Computational
  Statistics and Data Analysis}, In review.

  
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
   %%Further information on the package \pkg{spBayes} can be found at: \url{http://blue.fr.umn.edu/spatialBayes}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}
}

\examples{
\dontrun{
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
w <- rf.n200.dat$w

##############################
##Estimating all parameters
##############################
m.1 <- sp.lm(Y~1, coords=coords,
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100, sp.effects=TRUE)

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
##Fixing spatial and variance
##parameters
##############################
m.2 <- sp.lm(Y~1, coords=coords,
             fixed=list("phi"=0.6,"sigma.sq"=10, "tau.sq"=1),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100, sp.effects=TRUE)

print(summary(m.2$p.samples))
plot(m.2$p.samples)
 
##############################
##Modified predictive process
##############################
##Use some more observations
data(rf.n500.dat)

Y <- rf.n500.dat$Y
coords <- as.matrix(rf.n500.dat[,c("x.coords","y.coords")])
w <- rf.n500.dat$w

##############################
##Estimating all parameters
##############################
m.3 <- sp.lm(Y~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=2000, verbose=TRUE, n.report=100, sp.effects=TRUE)

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
legend(1.5,2.5, legend=c("Obs.", "Knots"), pch=c(1,19), cex=c(1,1), bg="white")
}
}

\keyword{misc}
