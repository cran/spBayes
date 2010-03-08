\name{mvLM}
\alias{mvLM}
\title{Function for fitting multivariate Bayesian regression models}

\description{
  The function \code{mvLM} fits a Gaussian multivariate Bayesian
  regression model.
}

\usage{
mvLM(formula, data = parent.frame(), starting, tuning, priors,
     n.samples, sub.samples, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{for a multivariate model with \eqn{q}{q} response variables, this is a list of univariate formulas.
    See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{mvLM} is called.  }

  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid list tags are \code{beta} and \code{L}.  The value portion of each tag
    is a vector of parameter's starting value. For \code{L}
    the vector is of length \eqn{\frac{(q^2-q)}{2}+q}{(q^2-q)/2+q}. Here,
    \code{L} holds the the lower-triangle elements in column major ordering of the Cholesky square root
    of the residual cross-covariance matrix. }
  \item{tuning}{a list with each tag corresponding to a
    parameter name. The only valid list tag is
    \code{L}. The value portion of each tag defines the step size of the proposal used in the
    Metropolis sampling. For \code{L}
    the vector is of length \eqn{\frac{(q^2-q)}{2}+q}{(q^2-q)/2+q}.  }
  \item{priors}{a list with each tag corresponding to a
    parameter name. The only valid list tag is \code{Psi.IW} (\code{Beta} priors are
    assumed \emph{flat}). \code{Psi} is assumed to follow an
    inverse-Wishart distribution. The
    hyperparameters of the inverse-Wishart are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{q\times q}{qxq} \emph{scale} matrix, respectively.}
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
  An object of class \code{mvLM}, which is a list with the following
  tags:
  \item{p.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis sampling
    acceptance rate. }
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\seealso{
  \code{\link{spMvLM}},\code{\link{spGGT}}
}

\examples{
\dontrun{
##Generate some data
set.seed(1)

n <- 100
q <- 3   
nltr <- q*(q-1)/2+q

L <- matrix(0,q,q)
L[lower.tri(L,TRUE)] <- rnorm(nltr,1,1)
Psi <- L\%*\%t(L)

X <- mkMvX(list(matrix(1,n,1), matrix(1,n,1), matrix(1,n,1)))
beta <- c(-1,0,1)

y <- mvrnorm(1, X\%*\%beta, diag(n)\%x\%Psi)

y.1 <- y[seq(1,length(y),q)]
y.2 <- y[seq(2,length(y),q)]
y.3 <- y[seq(3,length(y),q)]

n.samples <- 10000

##starting
L.starting <- diag(q)[lower.tri(diag(q),TRUE)]
L.tuning <- rep(0.03,nltr)

m.1 <- mvLM(list(y.1~1,y.2~1,y.3~1), starting=list("L"=L.starting), 
               tuning=list("L"=L.tuning), priors=list("Psi.IW"=list(q+1, diag(0.1,q))),
               n.samples=n.samples, verbose=TRUE, n.report=1000)

L.starting <- diag(5,q)[lower.tri(diag(q),TRUE)]
m.2 <- mvLM(list(y.1~1,y.2~1,y.3~1), starting=list("L"=L.starting), 
               tuning=list("L"=L.tuning), priors=list("Psi.IW"=list(q+1, diag(0.1,q))),
               n.samples=n.samples, verbose=TRUE, n.report=1000)

L.starting <- diag(10,q)[lower.tri(diag(q),TRUE)]
m.3 <- mvLM(list(y.1~1,y.2~1,y.3~1), starting=list("L"=L.starting), 
               tuning=list("L"=L.tuning), priors=list("Psi.IW"=list(q+1, diag(0.1,q))),
               n.samples=n.samples, verbose=TRUE, n.report=1000)

samps <- mcmc.list(m.1$p.samples, m.2$p.samples, m.3$p.samples)
gelman.plot(samps)

##Note that Psi[2,1], Psi[3,1], and Psi[3,2] are slow to converge.
##Try a bit of hand tuning or run the chains out longer.
n.samples <- 25000

##starting
L.starting <- diag(q)[lower.tri(diag(q),TRUE)]
L.tuning <- rep(0.03,nltr)

m.1 <- mvLM(list(y.1~1,y.2~1,y.3~1), starting=list("L"=L.starting), 
               tuning=list("L"=L.tuning),
               priors=list("Psi.IW"=list(q+1, diag(0.1,q))),
               n.samples=n.samples, verbose=TRUE, n.report=1000)

L.starting <- diag(5,q)[lower.tri(diag(q),TRUE)]
m.2 <- mvLM(list(y.1~1,y.2~1,y.3~1), starting=list("L"=L.starting), 
               tuning=list("L"=L.tuning),
               priors=list("Psi.IW"=list(q+1, diag(0.1,q))),
               n.samples=n.samples, verbose=TRUE, n.report=1000)

L.starting <- diag(10,q)[lower.tri(diag(q),TRUE)]
m.3 <- mvLM(list(y.1~1,y.2~1,y.3~1), starting=list("L"=L.starting), 
               tuning=list("L"=L.tuning),
               priors=list("Psi.IW"=list(q+1, diag(0.1,q))),
               n.samples=n.samples, verbose=TRUE, n.report=1000)

samps <- mcmc.list(m.1$p.samples, m.2$p.samples, m.3$p.samples)
gelman.plot(samps)

##parameters posterior summary
burn.in <- 10000
plot(window(samps, start=burn.in))

print(round(summary(window(samps, start=burn.in))$quantiles[,c(3,1,5)],2))
}
}

\keyword{misc}
