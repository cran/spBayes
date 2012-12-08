\name{bayesLMRef}
\alias{bayesLMRef}          
\title{Simple Bayesian linear model with non-informative priors }

\description{
  Given a \code{lm} object, the \code{bayesLMRef} function fits a
  simple Bayesian linear model with reference (non-informative) priors.  
}

\usage{
  bayesLMRef(lm.obj, n.samples, ...)
}

\arguments{
  \item{lm.obj}{an object returned by \code{lm}. }
  \item{n.samples}{the number of posterior samples to collect. }
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{bayesLMRef}, which is a list with at
  least the following tag:
  \item{p.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
}

\details{
  See page 355 in Gelman et al. (2004).
}

\references{
  Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004). Bayesian Data Analysis. 2nd ed. Boca Raton, FL: Chapman and Hall/CRC Press.

}
  
\author{
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\examples{
\dontrun{
set.seed(1)

n <- 100
X <- as.matrix(cbind(1, rnorm(n)))
B <- as.matrix(c(1,5))
tau.sq <- 0.1
y <- rnorm(n, X\%*\%B, sqrt(tau.sq))

lm.obj <- lm(y ~ X-1)

summary(lm.obj)

##Now with bayesLMRef
n.samples <- 500

m.1 <- bayesLMRef(lm.obj, n.samples)

summary(m.1$p.beta.tauSq.samples)
}
}

\keyword{model}