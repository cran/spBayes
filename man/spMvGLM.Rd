\name{spMvGLM}
\alias{spMvGLM}
\title{Function for fitting multivariate Bayesian generalized linear spatial regression models}

\description{
  The function \code{spMvGLM} fits multivariate stationary Bayesian
  generalized linear spatial regression models. Given a set of knots, \code{spMvGLM} fits a
  \emph{predictive process} model (see references below).
}

\usage{
spMvGLM(formula, family="binomial", data = parent.frame(), coords, knots,
      starting, tuning, priors, cov.model, n.samples, sub.samples,
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{for a multivariate model with \eqn{q}{q} response variables, this is a list of univariate formulas.
    See example below. }
  \item{family}{currently only supports \code{binomial} and
    \code{poisson} data using the logit and log link functions, respectively.}
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spMvGLM} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing)
    or a vector of length two or three with the first and second elements recording the
    number of columns and rows in the desired knot grid. The third,
    optional, element sets the offset of the outermost knots from the
    extent of the \code{coords} extent. } 
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{A}, \code{phi}, and \code{nu}.  The value portion of each tag
    is a vector of parameter's starting value. For \code{A} the vector
    is of length \eqn{\frac{q(q-q)}{2}+q}{q(q-q)/2+q}
    and \code{phi} and \code{nu} are of length \eqn{q}{q}. Here,
    \code{A} holds the the lower-triangle elements in column major ordering of the Cholesky square root
    of the spatial cross-covariance matrix. If the predictive
    process is used then \code{w} must be of length \eqn{qm}{qm}; otherwise,
    it must be of length \eqn{qn}{qn}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated. }
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{A},
    \code{phi}, \code{nu}, and \code{w}. The value portion of each tag defines the variance of the Metropolis normal proposal
    distribution. The tuning value for \code{beta} can be a
    vector of length \eqn{p}{p} or the lower-triangle of the
    \eqn{p\times p}{pxp} Cholesky square-root of the desired proposal variance matrix. For \code{A}, the vector is of length \eqn{\frac{q(q-q)}{2}+q}{q(q-q)/2+q}
    and \code{phi} and \code{nu} are of length \eqn{q}{q}. If the
    predictive process is used then \code{w} 
    must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which the value is repeated.  }
  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta.flat},
    \code{beta.normal}, \code{K.IW}, \code{Psi.IW}, \code{phi.unif}, and
    \code{nu.unif}. If \code{beta.normal} then covariate specific mean and variance hyperparameters are
    passed as the first and second list elements, respectively.  \code{K} is assumed to follow an
    inverse-Wishart distribution, whereas the spatial range \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Wishart are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{q\times q}{qxq} \emph{scale} matrix, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length \eqn{2\times q}{2xq} with consecutive elements representing the first
    and second elements corresponding to the lower and upper support in
  the order of the univariate models given in \code{formula}. }
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
  An object of class \code{spMvGLM}, which is a list with the following
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
  \code{\link{spGGT}},\code{\link{spMvLM}}, \code{\link{spMvGLM}}
}

\references{
    Finley, A.O., S. Banerjee, and R.E. McRoberts. (2008) A Bayesian
    approach to quantifying uncertainty in multi-source forest area
    estimates. \emph{Environmental and Ecological Statistics},
15:241--258.

  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian
  Predictive Process Models for Large Spatial Datasets. \emph{Journal of
  the Royal Statistical Society Series B}, 70:825--848.

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
################################
##Spatial multivariate poisson
################################
##Generate some data
n <- 100
q <- 3   
nltr <- q*(q-1)/2+q

coords <- cbind(runif(n,1,100),runif(n,1,100))

theta <- rep(3/50,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rnorm(nltr,1,1)
K <- A\%*\%t(A)

Psi <- diag(0,q)

c1 <- mvCovInvLogDet(coords=coords, knots=coords,
                     cov.model="exponential",
                     V=K, Psi=Psi, theta=theta,
                     modified.pp=TRUE, SWM=FALSE)

w <- mvrnorm(1,rep(0,nrow(c1$C)),c1$C)

X <- mkMvX(list(matrix(1,n,1), matrix(1,n,1), matrix(1,n,1)))
beta <- c(-1,0,1)
y <- rpois(n*q, exp(X\%*\%beta+w))

y.1 <- y[seq(1,length(y),q)]
y.2 <- y[seq(2,length(y),q)]
y.3 <- y[seq(3,length(y),q)]

##Specify starting values and collect samples. For
##a true analysis, several longer chains should be
##run.
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

beta.starting <- coefficients(glm(y~X-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~X-1, family="poisson"))))

n.samples <- 15000

m.1 <- spMvGLM(list(y.1~1,y.2~1,y.3~1), family="poisson", coords=coords,
               knots=c(8,8,0),
               starting=
               list("beta"=beta.starting, "phi"=rep(0.06,q),
                    "A"=A.starting, "w"=0), 
               tuning=
               list("beta"=beta.tuning, "phi"=rep(0.01,q),
                    "A"=rep(0.005,nltr), "w"=0.001),
               priors=
               list("beta.Flat", "phi.Unif"=rep(c(0.03, 0.3),q),
                    "K.IW"=list(q+1, diag(0.1,q))),
               cov.model="exponential",
               n.samples=n.samples, sub.sample=c(5000,n.samples,10),
               verbose=TRUE, n.report=500)

m.1$p.samples[,paste("phi_",1:q,sep="")] <-
  3/m.1$p.samples[,paste("phi_",1:q,sep="")]

print(summary(mcmc(m.1$p.samples)))

beta.hat <- colMeans(m.1$p.samples[,1:q])
w.hat <- rowMeans(m.1$sp.effects)

y.hat <- exp(X\%*\%beta.hat+w.hat)

y.hat.1 <- y.hat[seq(1,length(y.hat),q)]
y.hat.2 <- y.hat[seq(2,length(y.hat),q)]
y.hat.3 <- y.hat[seq(3,length(y.hat),q)]

##Take a look
par(mfrow=c(3,2))
surf <- mba.surf(cbind(coords,y.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.1, cex=1, col="blue")

surf <- mba.surf(cbind(coords,y.hat.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, add=TRUE)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat.1,0), cex=1, col="blue")

surf <- mba.surf(cbind(coords,y.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.2, cex=1, col="blue")

surf <- mba.surf(cbind(coords,y.hat.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat.2,0), cex=1, col="blue")

surf <- mba.surf(cbind(coords,y.3),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.3, cex=1, col="blue")

surf <- mba.surf(cbind(coords,y.hat.3),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat.3,0), cex=1, col="blue")
}
}

\keyword{misc}
