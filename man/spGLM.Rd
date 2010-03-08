\name{spGLM}
\alias{spGLM}
\title{Function for fitting univariate Bayesian generalized linear spatial regression models}

\description{
  The function \code{spGLM} fits univariate Bayesian
  generalized linear spatial regression models. Given a set of knots,
  \code{spGLM} will also fit a \emph{predictive process} model (see references below).
}

\usage{
spGLM(formula, family="binomial", weights, data = parent.frame(),
      coords, knots, amcmc, starting, tuning, priors, cov.model,
      modified.pp = TRUE, n.samples, sub.samples, verbose=TRUE,
      n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{family}{currently only supports \code{binomial} and
  \code{poisson} data using the logit and log link functions,
  respectively.}
  \item{weights}{an optional vector of weights to be used in the fitting
  process when \code{family} is \code{binomial}. Unless otherwise
  specified a vector of ones is assumed i.e., one trial at
  each location.}
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
    extent of the \code{coords}. }
  \item{amcmc}{a list with tags \code{n.batch}, \code{batch.length}, and
  \code{accept.rate}. Specifying this argument invokes an adaptive MCMC
  sampler (see Roberts and Rosenthal, 2006 and \code{adaptMetropGibbs} function).}
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{sigma.sq},
    \code{phi}, \code{nu}, and \code{w}.  The value portion of each tag is the parameter's starting value. If the predictive
    process is used then \code{w} must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated.}
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{sigma.sq},
    \code{phi}, \code{nu}, \code{w}, and \code{e}. The value portion of each tag defines the variance of the Metropolis normal proposal distribution. The tuning value for \code{beta} can be a
    vector of length \eqn{p}{p} or the lower-triangle of the
  \eqn{p\times p}{pxp} Cholesky square-root of the desired proposal variance matrix.
    If the predictive process is used then \code{w} must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which the value is repeated.  If a modified
  predictive process is not used then \code{e} is ignored.}
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
  \item{modified.pp}{a logical value indicating if the \emph{modified
      predictive process} should be used (see references below for
    details). Note, if a predictive process model is not used (i.e., \code{knots} is not specified) then
    this argument is ignored. }
  \item{n.samples}{the number of MCMC iterations. This is ignored if
  \code{amcmc} is specified. }
  \item{sub.samples}{a vector of length 3 that specifies \emph{start},
  \emph{end}, and \emph{thin}, respectively, of the MCMC samples. The
  default is \code{c(1, n.samples, 1)} (i.e., all samples). }
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report Metropolis acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\details{
  If a \code{binomial} model is specified the response vector is the
  number of successful trials at each location. Note, this is different
  from the \code{glm} argument which is passed as a rate. 
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
    acceptance rate. If \code{amcmc} is used then this will be a matrix of
    each parameter's acceptance rate at the end of each
    batch. Otherwise, the sampler is a Metropolis with a joint proposal
    of all parameters, and the acceptance rate is the average over all proposals.}
  \item{acceptance.w}{If this is a non-predictive process model and \code{amcmc} is used then this will be a matrix of
    each random spatial effects acceptance rate at the end of each
    batch. }
  \item{acceptance.w.str}{If this is a predictive process model and \code{amcmc} is used then this will be a matrix of
    each random spatial effects acceptance rate at the end of each batch. }
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

  Finley, A.O,. H. Sang, S. Banerjee, and A.E. Gelfand. (2008) Improving
  the performance of predictive process modeling for large datasets. \emph{Computational Statistics and Data Analysis}, DOI: 10.1016/j.csda.2008.09.008    

  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004) Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
  
  Roberts G.O. and Rosenthal J.S. (2006) Examples of Adaptive MCMC. \url{http://probability.ca/jeff/ftpdir/adaptex.pdf} Preprint.

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
set.seed(1)

n <- 100

coords <- cbind(runif(n,1,100),runif(n,1,100))

phi <- 3/50
sigma.sq <- 2

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- mvrnorm(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
y <- rpois(n, exp(x\%*\%beta+w))

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="poisson"))))

n.batch <- 500
batch.length <- 50
n.samples <- n.batch*batch.length

##Note tuning list is now optional

m.1 <- spGLM(y~1, family="poisson", coords=coords,
             starting=
             list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=
             list("beta"=0.1, "phi"=0.5, "sigma.sq"=0.1, "w"=0.1),
             priors=
             list("beta.Flat", "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=
             list("n.batch"=n.batch,"batch.length"=batch.length,
                  "accept.rate"=0.43),
             cov.model="exponential",
             n.samples=n.samples, sub.samples=c(1000,n.samples,10),
             verbose=TRUE, n.report=10)

##Just for fun check out the progression of the acceptance
##as it moves to 43\% (same can be seen for the random spatial effects).
plot(mcmc(t(m.1$acceptance)), density=FALSE, smooth=FALSE)

##Now parameter summaries, etc.
m.1$p.samples[,"phi"] <- 3/m.1$p.samples[,"phi"]

plot(mcmc(m.1$p.samples))
print(summary(mcmc(m.1$p.samples)))

beta.hat <- mean(m.1$p.samples[,"(Intercept)"])
w.hat <- rowMeans(m.1$sp.effects)

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
n <- 100

coords <- cbind(runif(n,1,100),runif(n,1,100))

phi <- 3/50
sigma.sq <- 2

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- mvrnorm(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
p <- 1/(1+exp(-(x\%*\%beta+w)))
y <- rbinom(n, 1, prob=p)

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="binomial"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="binomial"))))
            
n.batch <- 500
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y~1, family="binomial", coords=coords, 
             starting=
             list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=
             list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.1, "w"=0.01),
             priors=
             list("beta.Normal"=list(0,10), "phi.Unif"=c(0.03, 0.3),
                  "sigma.sq.IG"=c(2, 1)),
             amcmc=
             list("n.batch"=n.batch,"batch.length"=batch.length,
                  "accept.rate"=0.43),
             cov.model="exponential",
             n.samples=n.samples, sub.samples=c(1000,n.samples,10),
             verbose=TRUE, n.report=10)

m.1$p.samples[,"phi"] <- 3/m.1$p.samples[,"phi"]

print(summary(mcmc(m.1$p.samples)))

beta.hat <- mean(m.1$p.samples[,"(Intercept)"])
w.hat <- rowMeans(m.1$sp.effects)

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
