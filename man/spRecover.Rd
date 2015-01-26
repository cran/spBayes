\name{spRecover}
\alias{spRecover}
\title{Function for recovering regression coefficients and spatial
  random effects for \code{\link{spLM}}, \code{\link{spMvLM}}, and
  \code{\link{spMisalignLM}} using composition sampling
}

\description{
 Function for recovering regression coefficients and spatial random
 effects for \code{\link{spLM}}, \code{\link{spMvLM}}, and
  \code{\link{spMisalignLM}} using composition sampling.
}

\usage{
spRecover(sp.obj, get.beta=TRUE, get.w=TRUE, start=1, end, thin=1,
          verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{spLM}},
  \code{\link{spMvLM}}, or \code{\link{spMisalignLM}}.}
  \item{get.beta}{if \code{TRUE}, regression coefficients will be recovered.}
  \item{get.w}{if \code{TRUE}, spatial random effects will be recovered.}
  \item{start}{specifies the first sample included in the composition sampling.}
  \item{end}{specifies the last sample included in the composition.
             The default is to use all posterior samples in \code{sp.obj}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
              samples between \code{start} and \code{end}.  For example, if \code{thin = 10} then 1 in 10 samples are considered between \code{start} and
              \code{end}. }
      \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report sampling progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  The input \code{sp.obj} with posterior samples of regression coefficients and/or spatial random effects appended. 
  tags:
    \item{p.theta.recover.samples}{those \code{p.theta.samples} used in the composition sampling. }
    \item{p.beta.recover.samples}{a \code{coda} object of regression coefficients posterior samples. }
    \item{p.w.recover.samples}{a \code{coda} object of spatial random
      effects posterior samples. Rows correspond to locations' random
      effects and columns are posterior samples.}  
}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical
  modeling and analysis for spatial data. Chapman and Hall/CRC Press,
  Boca Raton, FL.

    Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 50
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)
sigma.sq <- 10
tau.sq <- 0.01
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X\%*\%B + w, sqrt(tau.sq))

n.samples <- 1000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
               "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 0.01))
cov.model <- "exponential"

m.1 <- spLM(y~X-1, coords=coords, starting=starting, tuning=tuning,
            priors=priors, cov.model=cov.model, n.samples=n.samples)

m.1 <- spRecover(m.1, start=0.5*n.samples, thin=2)

summary(window(m.1$p.beta.recover.samples))

w.hat <- apply(m.1$p.w.recover.samples, 1, mean)
plot(w, w.hat, xlab="Observed w", ylab="Fitted w")
}
}

\keyword{utilities}