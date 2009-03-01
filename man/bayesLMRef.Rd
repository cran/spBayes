\name{bayesLMRef}
\alias{bayesLMRef}          
\title{Simple Bayesian linear model with non-informative priors }

\description{
  Given an \code{lm} object, the \code{bayesLMRef} function fits a
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
data(FORMGMT.dat)

lm.obj <- lm(Y ~ X1+X2+X3+X4+X5+X6, data = FORMGMT.dat)

summary(lm.obj)

##Now with bayesLMRef
n.samples <- 500

m.1 <- bayesLMRef(lm.obj, n.samples)

summary(m.1$p.sample)

}

\keyword{misc}
