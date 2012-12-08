\name{spDiag}
\alias{spDiag}          
\title{Model fit diagnostics}

\description{
  The function \code{spDiag} calculates DIC, GP, GRS, and associated
  statistics given a \code{\link{spLM}}, \code{\link{spMvLM}},
  \code{\link{spGLM}}, or \code{\link{spMvGLM}} object.
}

\usage{
  spDiag(sp.obj, start=1, end, thin=1, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{sp.obj}{an object returned by  \code{\link{spLM}},
  \code{\link{spMvLM}}, \code{\link{spGLM}}, or \code{\link{spMvGLM}}.}
  \item{start}{specifies the first sample included in the computation. The \code{start}, \code{end}, and \code{thin} arguments only apply to \code{\link{spGLM}} or
  \code{\link{spMvGLM}} objects. Sub-sampling for \code{\link{spLM}} and \code{\link{spMvLM}} is controlled using \code{\link{spRecover}} which
  must be called prior to \code{spDiag}.}
  \item{end}{specifies the last sample included  in the computation.
    The default is to use all posterior samples in \code{sp.obj}. See
  \code{start} argument description.}
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10} then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{verbose}{if \code{TRUE} calculation progress is printed to the
    screen; otherwise, nothing is printed to the screen.  }
  \item{n.report}{the interval to report progress. } 
  \item{...}{currently no additional arguments.  }
}

\value{
  A list with some of the following tags:
  \item{DIC}{a matrix holding DIC and associated statistics, see
  Banerjee et al. (2004) for details.}
  \item{GP}{a matrix holding GP and associated statistics, see Gelfand
    and Ghosh (1998) for details. }
  \item{GRS}{a scoring rule, see Equation 27 in Gneiting and Raftery
  (2007) for details.}
}


\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical
  modeling and analysis for spatial data. Chapman and Hall/CRC Press,
  Boca Raton,Fla.
  
  Gelfand A.E. and Ghosh, S.K. (1998). Model choice: a minimum posterior
  predictive loss approach. \emph{Biometrika}. 85:1-11.

  Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring rules, prediction, and estimation.
  \emph{Journal of the American Statistical Association}. 102:359-378. 
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}.
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

n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- 2
tau.sq <- 0.1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X\%*\%B + w, sqrt(tau.sq))

n.samples <- 1000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

##too restrictive of prior on beta
priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1,p)),
                 "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

##more reasonable prior for beta
priors.2 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

m.1 <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

m.2 <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.2, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

##non-spatial model
m.3 <- spLM(y~X-1, n.samples=n.samples, verbose=verbose, n.report=n.report)

burn.in <- 0.5*n.samples

##recover beta and spatial random effects
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)
m.2 <- spRecover(m.2, start=burn.in, verbose=FALSE)

##lower is better for DIC, GPD, and GRS
print(spDiag(m.1))
print(spDiag(m.2))
print(spDiag(m.3))
}
}

\keyword{utilities}