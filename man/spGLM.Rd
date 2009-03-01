\name{spGLM}
\alias{spGLM}
\title{Function for fitting univariate Bayesian generalized linear spatial regression models}

\description{
  The function \code{spGLM} fits univariate stationary Bayesian
  generalized linear spatial regression models. Given a set of knots, \code{spGLM} fits a
  \emph{predictive process} model (see references below).
}

\usage{
spGLM(formula, family="binomial", data = parent.frame(), coords, knots,
      starting, tuning, priors, cov.model,
      n.samples, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{family}{currently only supports \code{binomial} and
  \code{poisson} data using the logit and log link functions, respectively.}
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spGLM} is called.  }
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
    \code{phi}, \code{nu}, and \code{w}.  The value portion of each tag is the parameter's starting value. If the predictive
    process is used then \code{w} must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated.}
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{sigma.sq},
    \code{phi}, \code{nu}, and \code{w}. The value portion of each tag defines the variance of the Metropolis normal proposal distribution. The tuning value for \code{beta} can be a
    vector of length \eqn{p}{p} or the lower-triangle of the
  \eqn{p\times p}{pxp} Cholesky square-root of the desired proposal variance matrix.
    If the predictive process is used then \code{w} must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which the value is repeated.  }
  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta.flat}, \code{beta.normal}, \code{sigma.sq.ig},
    \code{phi.unif}, and \code{nu.unif}. \code{simga.sq} is assumed to follow an
    inverse-Gamma distribution, whereas the spatial range \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform
    distributions. If \code{beta.normal} then covariate specific mean and variance hyperparameters are
    passed as the first and second list elements, respectively. The
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
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report Metropolis acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{spGLM}, which is a list with the following
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
  \code{\link{spGGT}}, \code{\link{spMvLM}}, \code{\link{spMvGLM}}
}

\references{
  Finley, A.O., S. Banerjee, and R.E. McRoberts. (2008) A Bayesian approach to quantifying uncertainty in multi-source forest area estimates. \emph{Environmental and Ecological Statistics}, 15:241--258.
  
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.

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
###########################
##Spatial poisson
###########################

##Generate count data
n <- 100

coords <- cbind(runif(n,1,100),runif(n,1,100))

phi <- 3/75
sigma.sq <- 2

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- mvrnorm(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
y <- rpois(n, exp(x\%*\%beta+w))

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="poisson"))))
            
n.samples <- 25000

m.1 <- spGLM(y~1, family="poisson", coords=coords, 
             starting=
             list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=
             list("beta"=0.1, "phi"=0.01, "sigma.sq"=0.01, "w"=0.005),
             priors=
             list("beta.Flat", "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=n.samples, verbose=TRUE, n.report=500)

m.1$p.samples[,"phi"] <- 3/m.1$p.samples[,"phi"]


burn.in <- 0.75*n.samples
print(summary(mcmc(m.1$p.samples[burn.in:n.samples,])))

beta.hat <- mean(m.1$p.samples[burn.in:n.samples,1])
w.hat <- rowMeans(m.1$sp.effects[,burn.in:n.samples])

y.hat <-exp(x\%*\%beta.hat+w.hat)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="obs")
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y, cex=1, col="blue")

surf <- mba.surf(cbind(coords,y.hat),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat,0), cex=1, col="blue")

###########################
##Spatial logistic
###########################

##Generate binary data
n <- 50

coords <- cbind(runif(n,1,100),runif(n,1,100))

phi <- 3/50
sigma.sq <- 1

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- mvrnorm(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
p <- 1/(1+exp(-(x\%*\%beta+w)))
y <- rbinom(n, 1, prob=p)

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="binomial"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="binomial"))))
            
n.samples <- 50000

m.1 <- spGLM(y~1, family="binomial", coords=coords, 
             starting=
             list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=
             list("beta"=beta.tuning, "phi"=0.1, "sigma.sq"=0.1, "w"=0.01),
             priors=
             list("beta.Normal"=list(0,10), "phi.Unif"=c(0.03, 0.3),
                  "sigma.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=n.samples, verbose=TRUE, n.report=500)

m.1$p.samples[,"phi"] <- 3/m.1$p.samples[,"phi"]

burn.in <- 0.75*n.samples
print(summary(mcmc(m.1$p.samples[burn.in:n.samples,])))

beta.hat <- mean(m.1$p.samples[burn.in:n.samples,1])
w.hat <- rowMeans(m.1$sp.effects[,burn.in:n.samples])

y.hat <- 1/(1+exp(-(x\%*\%beta.hat+w.hat)))

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed")
contour(surf, add=TRUE)
points(coords[y==1,], pch=19, cex=1)
points(coords[y==0,], cex=1)

surf <- mba.surf(cbind(coords,y.hat),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted probabilities")
contour(surf, add=TRUE)
points(coords[y==1,], pch=19, cex=1)
points(coords[y==0,], cex=1)

}
}

\keyword{misc}