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
  \item{p.samples}{a \code{coda} object of posterior samples for the parameters. }
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
w <-  mvrnorm(1, rep(0, n), s)
Y <- mvrnorm(1, rep(mu, n) + w, tausq*diag(n))
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
chain.1$p.samples[,2] <- exp(chain.1$p.samples[,2])
chain.1$p.samples[,3] <- exp(chain.1$p.samples[,3])
chain.1$p.samples[,4] <- 3/logit.inv(chain.1$p.samples[,4], a.phi, b.phi)

chain.2$p.samples[,2] <- exp(chain.2$p.samples[,2])
chain.2$p.samples[,3] <- exp(chain.2$p.samples[,3])
chain.2$p.samples[,4] <- 3/logit.inv(chain.2$p.samples[,4], a.phi, b.phi)

par.names <- c("mu", "sigma.sq", "tau.sq", "effective range (-log(0.05)/phi)")
colnames(chain.1$p.samples) <- par.names
colnames(chain.2$p.samples) <- par.names

##Discard burn.in and plot and do some convergence diagnostics
chains <- mcmc.list(mcmc(chain.1$p.samples), mcmc(chain.2$p.samples))
plot(window(chains, start=as.integer(0.5*n.batch*batch.length)))

gelman.diag(chains)

##########################
##Example of fitting a
##a non-linear model
##########################
set.seed(1)

##Simulate some data
par <- c(0.1,0.1,3,1)
tau.sq <- 0.01

fn <- function(par,x){
  par[1]+par[2]^(1-(log(x/par[3])/par[4])^2)
}

n <- 200
x <- seq(1,10,0.1)
y <- rnorm(length(x), fn(par,x), sqrt(tau.sq))

##Define the log target density
ltd <- function(theta, x, y, fn, IG.a, IG.b){

  theta[2:5] <- exp(theta[2:5])
  tau.sq <- theta[5]

  y.hat <- fn(theta[1:4],x)

  ##likelihood
  logl <- -(n/2)*log(tau.sq)-(1/(2*tau.sq))*sum((y-y.hat)^2)

  ##priors IG on tau.sq and the rest flat
  logl <- logl-(IG.a+1)*log(tau.sq)-IG.b/tau.sq

  ##Jacobian adjustments
  logl <- logl+sum(log(theta[2:5]))
  
  return(logl)  
}

##Give some reasonable starting values,
##note the transformation
starting <- c(0.5, log(0.5), log(5), log(3), log(1))

m.1 <- adaptMetropGibbs(ltd, starting=starting,
                        batch=1200, batch.length=25, 
                        x=x, y=y, fn=fn, IG.a=2, IG.b=0.01)

##Back transform
m.1$p.samples[,2:5] <- exp(m.1$p.samples[,2:5])

##Summary with 95% prediction interval
burn.in <- 10000

fn.pred <- function(par,x){
  rnorm(length(x), fn(par[1:4],x), sqrt(par[5]))
}

post.curves <-
  t(apply(m.1$p.samples[burn.in:nrow(m.1$p.samples),], 1, fn.pred, x))

post.curves.quants <- summary(mcmc(post.curves))$quantiles

plot(x, y, pch=19, xlab="x", ylab="f(x)")
lines(x, post.curves.quants[,1], lty="dashed", col="blue")
lines(x, post.curves.quants[,3])
lines(x, post.curves.quants[,5], lty="dashed", col="blue")

########################################
##Now some real data. There are way too
##few data for this model so we use a
##slightly informative normal prior on
##the parameters.
########################################

fn <- function(par,x){
  par[1]-(par[2]^(1-(log(x/par[3])/par[4])^2))
}

x <- c(5,6,7,17,17,18,39,55,60)
y <- c(47.54,62.10,37.29,24.74,22.64,32.60,11.16,9.65,14.83)
n <- length(x)

##Define the log target density
ltd <- function(theta, x, y, fn, IG.a, IG.b){

  theta[2:5] <- exp(theta[2:5])
  tau.sq <- theta[5]

  y.hat <- fn(theta[1:4],x)

  ##likelihood
  logl <- -(n/2)*log(tau.sq)-(1/(2*tau.sq))*sum((y-y.hat)^2)

  ##priors IG on tau.sq and the rest normal
  logl <- logl-(IG.a+1)*log(tau.sq)-IG.b/tau.sq+
    sum(dnorm(theta[1:4], rep(0,4), 500, log = TRUE))
  
  ##Jacobian adjustments
  logl <- logl+sum(log(theta[2:5]))
  
  return(logl)  
}

##Give some reasonable starting values,
##note the transformation
starting <- c(50, log(40), log(40), log(2), log(0.1))

m.1 <- adaptMetropGibbs(ltd, starting=starting,
                        batch=1500, batch.length=25, 
                        x=x, y=y, fn=fn, IG.a=2, IG.b=10)

##Back transform
m.1$p.samples[,2:5] <- exp(m.1$p.samples[,2:5])

##Summary with 95% prediction interval
burn.in <- 20000

x.pred <- seq(0, 60, seq=0.1)

fn.pred <- function(par,x){
  rnorm(length(x), fn(par[1:4],x), sqrt(par[5]))
}

post.curves <-
  t(apply(m.1$p.samples[burn.in:nrow(m.1$p.samples),], 1, fn.pred, x.pred))

post.curves.quants <- summary(mcmc(post.curves))$quantiles

plot(x, y, pch=19, xlab="x", ylab="f(x)", xlim=c(0,60), ylim=c(0,100))
lines(x.pred, post.curves.quants[,1], lty="dashed", col="blue")
lines(x.pred, post.curves.quants[,3])
lines(x.pred, post.curves.quants[,5], lty="dashed", col="blue")

}
}

\keyword{misc}
