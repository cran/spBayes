\name{spMvLM}
\alias{spMvLM}
\title{Function for fitting multivariate Bayesian spatial regression models}

\description{
  The function \code{spMvLM} fits Gaussian multivariate Bayesian
  spatial regression models. Given a set of knots, \code{spMvLM} will
  also fit a \emph{predictive process} model (see references below).
}

\usage{
spMvLM(formula, data = parent.frame(), coords, knots,
      starting, sp.tuning, priors, cov.model,
      modified.pp = TRUE, n.samples, sub.samples,
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{for a multivariate model with \eqn{q}{q} response variables, this is a list of univariate formulas.
    See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spMvLM} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing)
    or a vector of length two or three with the first and second elements recording the
    number of columns and rows in the desired knot grid. The third,
    optional, element sets the offset of the outermost knots from the
    extent of the \code{coords} extent. } 
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta}, \code{A},
    \code{L}, \code{phi}, and \code{nu}.  The value portion of each tag
    is a vector of parameter's starting value. For \code{A} and \code{L}
    the vectors are of length \eqn{\frac{(q^2-q)}{2}+q}{(q^2-q)/2+q}
    and \code{phi} and \code{nu} are of length \eqn{q}{q}. Here,
    \code{A} and \code{L} hold the the lower-triangle elements in column major ordering of the Cholesky square root
    of the spatial and non-spatial cross-covariance matrices, respectively. }
  \item{sp.tuning}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{A},
    \code{L}, \code{phi}, and \code{nu}. The value portion of each
    of each tag defines the step size of the proposal used in the
    Metropolis sampling. For \code{A} and \code{L}
    the vectors are of length  \eqn{\frac{(q^2-q)}{2}+q}{(q^2-q)/2+q}
    and \code{phi} and \code{nu} are of length \eqn{q}{q}.  }
  \item{modified.pp}{a logical value indicating if the \emph{modified
    predictive process} should be used (see references below for
    details). Note, if a predictive process model is not used (i.e., \code{knots} is not specified) then
    this argument is ignored. }
  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{K.IW}, \code{Psi.IW}, \code{phi.unif}, and \code{nu.unif} (\code{Beta} priors are
    assumed \emph{flat}). Variance parameters, \code{K} and
    \code{Psi}, are assumed to follow an
    inverse-Wishart distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Wishart are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{q\times q}{qxq} \emph{scale} matrix, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length \eqn{2\times
  q}{2xq} with consecutive elements representing the first
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
  An object of class \code{spMvLM}, which is a list with the following
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
  \code{\link{spMvLM}},\code{\link{spGGT}}
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
###################################
##Multivariate spatial regression
###################################

##Generate some data
n <- 50 ##observed
q <- 3   
nltr <- q*(q-1)/2+q

coords <- cbind(runif(n),runif(n))
theta <- rep(3/0.5,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rnorm(nltr, 5, 1)
K <- A\%*\%t(A)

Psi <- diag(1,q)

c1 <- mvCovInvLogDet(coords=coords, knots=coords,
                     cov.model="exponential",
                     V=K, Psi=Psi, theta=theta,
                     modified.pp=FALSE, SWM=FALSE)

w <- mvrnorm(1,rep(0,nrow(c1$C)),c1$C)

w.1 <- w[seq(1,length(w),q)]
w.2 <- w[seq(2,length(w),q)]
w.3 <- w[seq(3,length(w),q)]

##Specify starting values and collect samples
K.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
Psi.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

n.samples <- 5000

m.1 <- spMvLM(list(w.1~1,w.2~1,w.3~1), coords=coords,
              starting=list("beta"=rep(1,q), "phi"=rep(3/0.5,q),
                "A"=K.starting, "L"=Psi.starting), 
              sp.tuning=list("phi"=rep(0.1,q),
                "A"=rep(0.1,nltr), "L"=rep(0.1,nltr)),
              priors=list("phi.Unif"=rep(c(3/1,3/0.1),q),
                "K.IW"=list(q+1, diag(1,q)), "Psi.IW"=list(q+1, diag(1,q))),
              modified.pp=TRUE, cov.model="exponential",
              n.samples=n.samples, sub.samples=c(2500,n.samples,5),
              verbose=TRUE, n.report=100)

m.1$p.samples[,paste("phi_",1:q,sep="")] <-
  3/m.1$p.samples[,paste("phi_",1:q,sep="")]

print(summary(mcmc(m.1$p.samples)))

w.hat <- rowMeans(m.1$sp.effects)
w.hat.1 <- w.hat[seq(1,length(w.hat),q)]
w.hat.2 <- w.hat[seq(2,length(w.hat),q)]
w.hat.3 <- w.hat[seq(3,length(w.hat),q)]

##Take a look
par(mfrow=c(3,2))
surf <- mba.surf(cbind(coords,w.1),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Observed 1"); contour(surf, add=TRUE)

surf <- mba.surf(cbind(coords,w.hat.1),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Fitted 1"); contour(surf, add=TRUE)
points(m.1$knot.coords, pch=19, cex=1)

surf <- mba.surf(cbind(coords,w.2),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Observed 2"); contour(surf, add=TRUE)

surf <- mba.surf(cbind(coords,w.hat.2),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Fitted 2"); contour(surf, add=TRUE)
points(m.1$knot.coords, pch=19, cex=1)

surf <- mba.surf(cbind(coords,w.3),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Observed 3"); contour(surf, add=TRUE)

surf <- mba.surf(cbind(coords,w.hat.3),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Fitted 3"); contour(surf, add=TRUE)
points(m.1$knot.coords, pch=19, cex=1)
}
}

\keyword{misc}
