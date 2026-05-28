\name{adaptMetropGibbs}
\alias{adaptMetropGibbs}
\title{Adaptive Metropolis within Gibbs algorithm}

\description{
  Markov chain Monte Carlo for continuous random vector using an
  adaptive Metropolis within Gibbs algorithm.
}

\usage{
adaptMetropGibbs(ltd, starting, tuning=1, accept.rate=0.44,
                 batch = 1, batch.length=25, report=100,
                 verbose=TRUE, ...)
}

\arguments{
  \item{ltd}{an \R function that evaluates the log target density of the
    desired equilibrium distribution of the Markov chain. First argument is the
    starting value vector of the Markov chain. Pass variables used in
    the \code{ltd} via the \ldots argument of \code{aMetropGibbs}.}
  \item{starting}{a real vector of parameter starting values. }
  \item{tuning}{ a scalar or vector of initial Metropolis tuning values. The vector must be of \code{length(starting)}. If a scalar is passed then it is expanded to \code{length(starting)}.}
  \item{accept.rate}{ a scalar or vector of target Metropolis acceptance
  rates. The vector must be of \code{length(starting)}. If a scalar is passed then it is expanded to \code{length(starting)}.}
  \item{batch}{the number of batches. } 
  \item{batch.length}{the number of sampler iterations in each batch. }
  \item{report}{the number of batches between acceptance rate reports. }
  \item{verbose}{if \code{TRUE}, progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{...}{currently no additional arguments.  }
}

\value{
  A list with the following tags:
  \item{p.theta.samples}{a \code{coda} object of posterior samples for the parameters. }
  \item{acceptance}{the Metropolis acceptance rate at the end of each
    batch. }
  \item{ltd}{\code{ltd} }
  \item{accept.rate}{\code{accept.rate} }
  \item{batch}{\code{batch} }
  \item{batch.length}{\code{batch.length} }
  \item{proc.time}{the elapsed CPU and wall time (in seconds). }
}

\references{
Roberts G.O. and Rosenthal J.S. (2006). Examples of Adaptive MCMC. \url{http://probability.ca/jeff/ftpdir/adaptex.pdf} Preprint.

Rosenthal J.S. (2007). AMCMC: An R interface for adaptive MCMC. \emph{Computational
  Statistics and Data Analysis}. 51:5467-5470.
}

\note{This function is a rework of Rosenthal (2007) with some added
  niceties. }

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}
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

###########################
##Fit a spatial regression
###########################
set.seed(1)
n <- 50
x <- runif(n, 0, 100)
y <- runif(n, 0, 100)

D <- as.matrix(dist(cbind(x, y)))

phi <- 3/50
sigmasq <- 50
tausq <- 20
mu <- 150

s <- (sigmasq*exp(-phi*D))
w <-  rmvn(1, rep(0, n), s)
Y <- rmvn(1, rep(mu, n) + w, tausq*diag(n))
X <- as.matrix(rep(1, length(Y)))

##Priors
##IG sigma^2 and tau^2
a.sig <- 2 
b.sig <- 100
a.tau <- 2
b.tau <- 100

##Unif phi
a.phi <- 3/100
b.phi <- 3/1

##Functions used to transform phi to continuous support.
logit <- function(theta, a, b){log((theta-a)/(b-theta))}
logit.inv <- function(z, a, b){b-(b-a)/(1+exp(z))}

##Metrop. target
target <- function(theta){
  
  mu.cand <- theta[1]
  sigmasq.cand <- exp(theta[2])
  tausq.cand <- exp(theta[3])
  phi.cand <- logit.inv(theta[4], a.phi, b.phi)

  Sigma <- sigmasq.cand*exp(-phi.cand*D)+tausq.cand*diag(n)
  SigmaInv <- chol2inv(chol(Sigma))
  logDetSigma <- determinant(Sigma, log=TRUE)$modulus[1]
  
  out <- (
          ##Priors
          -(a.sig+1)*log(sigmasq.cand) - b.sig/sigmasq.cand
          -(a.tau+1)*log(tausq.cand) - b.tau/tausq.cand
          ##Jacobians
          +log(sigmasq.cand) + log(tausq.cand) 
          +log(phi.cand - a.phi) + log(b.phi -phi.cand) 
          ##Likelihood
          -0.5*logDetSigma-0.5*(t(Y-X\%*\%mu.cand)\%*\%SigmaInv\%*\%(Y-X\%*\%mu.cand))
          )
  
  return(out)
}


##Run a couple chains
n.batch <- 500
batch.length <- 25

inits <- c(0, log(1), log(1), logit(3/10, a.phi, b.phi))
chain.1 <- adaptMetropGibbs(ltd=target, starting=inits,
                            batch=n.batch, batch.length=batch.length, report=100)

inits <- c(500, log(100), log(100), logit(3/90, a.phi, b.phi))
chain.2 <- adaptMetropGibbs(ltd=target, starting=inits,
                            batch=n.batch, batch.length=batch.length, report=100)

##Check out acceptance rate just for fun
plot(mcmc.list(mcmc(chain.1$acceptance), mcmc(chain.2$acceptance)))

##Back transform
chain.1$p.theta.samples[,2] <- exp(chain.1$p.theta.samples[,2])
chain.1$p.theta.samples[,3] <- exp(chain.1$p.theta.samples[,3])
chain.1$p.theta.samples[,4] <- 3/logit.inv(chain.1$p.theta.samples[,4], a.phi, b.phi)

chain.2$p.theta.samples[,2] <- exp(chain.2$p.theta.samples[,2])
chain.2$p.theta.samples[,3] <- exp(chain.2$p.theta.samples[,3])
chain.2$p.theta.samples[,4] <- 3/logit.inv(chain.2$p.theta.samples[,4], a.phi, b.phi)

par.names <- c("mu", "sigma.sq", "tau.sq", "effective range (-log(0.05)/phi)")
colnames(chain.1$p.theta.samples) <- par.names
colnames(chain.2$p.theta.samples) <- par.names

##Discard burn.in and plot and do some convergence diagnostics
chains <- mcmc.list(mcmc(chain.1$p.theta.samples), mcmc(chain.2$p.theta.samples))
plot(window(chains, start=as.integer(0.5*n.batch*batch.length)))

gelman.diag(chains)

##########################
##Example of fitting a
##a non-linear model
##########################
##Example of fitting a non-linear model
set.seed(1)

########################################################
##Simulate some data.
########################################################
a <- 0.1 #-Inf < a < Inf
b <- 0.1 #b > 0
c <- 0.2 #c > 0
tau.sq <- 0.1 #tau.sq > 0

fn <- function(a,b,c,x){
  a+b*exp(x/c)
}

n <- 200
x <- seq(0,1,0.01)
y <- rnorm(length(x), fn(a,b,c,x), sqrt(tau.sq))

##check out your data
plot(x, y)

########################################################
##The log target density
########################################################
##Define the log target density used in the Metrop.
ltd <- function(theta){

  ##extract and transform as needed
  a <- theta[1]
  b <- exp(theta[2])
  c <- exp(theta[3])
  tau.sq <- exp(theta[4])

  y.hat <- fn(a, b, c, x)

  ##likelihood
  logl <- sum(dnorm(y, y.hat, sqrt(tau.sq), log=TRUE))

  ##priors IG on tau.sq and normal on a and transformed b, c, d
  logl <- (logl
           -(IG.a+1)*log(tau.sq)-IG.b/tau.sq
           +sum(dnorm(theta[1:3], N.mu, N.v, log=TRUE))
           )
  
  ##Jacobian adjustment for tau.sq
  logl <- logl+log(tau.sq)
  
  return(logl)  
}

########################################################
##The rest
########################################################

##Priors
IG.a <- 2
IG.b <- 0.01

N.mu <- 0
N.v <- 10

theta.tuning <- c(0.01, 0.01, 0.005, 0.01)

##Run three chains with different starting values
n.batch <- 1000
batch.length <- 25

theta.starting <- c(0, log(0.01), log(0.6), log(0.01))
run.1 <- adaptMetropGibbs(ltd=ltd, starting=theta.starting, tuning=theta.tuning,
                          batch=n.batch, batch.length=batch.length, report=100)

theta.starting <- c(1.5, log(0.05), log(0.5), log(0.05))
run.2 <- adaptMetropGibbs(ltd=ltd, starting=theta.starting, tuning=theta.tuning,
                          batch=n.batch, batch.length=batch.length, report=100)

theta.starting <- c(-1.5, log(0.1), log(0.4), log(0.1))
run.3 <- adaptMetropGibbs(ltd=ltd, starting=theta.starting, tuning=theta.tuning,
                          batch=n.batch, batch.length=batch.length, report=100)

##Back transform
samples.1 <- cbind(run.1$p.theta.samples[,1], exp(run.1$p.theta.samples[,2:4]))
samples.2 <- cbind(run.2$p.theta.samples[,1], exp(run.2$p.theta.samples[,2:4]))
samples.3 <- cbind(run.3$p.theta.samples[,1], exp(run.3$p.theta.samples[,2:4]))

samples <- mcmc.list(mcmc(samples.1), mcmc(samples.2), mcmc(samples.3))

##Summary 
plot(samples, density=FALSE)
gelman.plot(samples)

burn.in <- 5000

fn.pred <- function(theta,x){
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  tau.sq <- theta[4]
  
  rnorm(length(x), fn(a,b,c,x), sqrt(tau.sq))
}

post.curves <- t(apply(samples.1[burn.in:nrow(samples.1),], 1, fn.pred, x))

post.curves.quants <- summary(mcmc(post.curves))$quantiles

plot(x, y, pch=19, xlab="x", ylab="f(x)")
lines(x, post.curves.quants[,1], lty="dashed", col="blue")
lines(x, post.curves.quants[,3])
lines(x, post.curves.quants[,5], lty="dashed", col="blue")


}
}
\keyword{model}