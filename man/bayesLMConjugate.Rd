\name{bayesLMConjugate}
\alias{bayesLMConjugate}          
\title{Simple Bayesian linear model via the Normal/inverse-Gamma conjugate }

\description{
  Given an \code{lm} object, the \code{bayesLMConjugate} function fits a
  simple Bayesian linear model with Normal and inverse-Gamma priors.  
}

\usage{
  bayesLMConjugate(formula, data = parent.frame(), n.samples,
                     beta.prior.mean, beta.prior.precision,
                     prior.shape, prior.rate, ...)
}

\arguments{
  \item{formula}{for a univariate model, this is a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spLM} is called.  }
  \item{n.samples}{the number of posterior samples to collect. }
  \item{beta.prior.mean}{\eqn{\beta}{beta} multivariate normal mean vector hyperprior. }
  \item{beta.prior.precision}{\eqn{\beta}{beta} multivariate normal precision
    matrix hyperprior. }
  \item{prior.shape}{\eqn{\sigma^2}{sigma.sq} inverse-Gamma shape
    hyperprior. }
  \item{prior.rate}{\eqn{\sigma^2}{sigma.sq} inverse-Gamma 1/scale
    hyperprior. }
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{bayesLMConjugate}, which is a list with at
  least the following tag:
  \item{p.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
}

%\references{
%  \url{http://blue.fr.umn.edu/spatialBayes}
%}
  
\author{
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\examples{
\dontrun{

data(FORMGMT.dat)

n = nrow(FORMGMT.dat)
p = 7 ##an intercept and six covariates

n.samples <- 500

## Below we demonstrate the conjugate function in the special case
## with improper priors. The results are the same as for the above,
## up to MC error. 
beta.prior.mean <- rep(0, times=p)
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

prior.shape <- -p/2
prior.rate <- 0

m.1 <-
  bayesLMConjugate(Y ~ X1+X2+X3+X4+X5+X6, data = FORMGMT.dat,
                     n.samples, beta.prior.mean,
                     beta.prior.precision,
                     prior.shape, prior.rate)

summary(m.1$p.samples)
}
}

\keyword{misc}
