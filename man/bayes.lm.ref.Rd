\name{bayes.lm.ref}
\alias{bayes.lm.ref}          
\title{Simple Bayesian linear model with non-informative priors }

\description{
  Given an \code{lm} object, the \code{bayes.lm.ref} function fits a
  simple Bayesian linear model with reference (non-informative) priors.  
}

\usage{
  bayes.lm.ref(lm.obj, n.samples, ...)
}

\arguments{
  \item{lm.obj}{an object returned by \code{lm}. }
  \item{n.samples}{the number of posterior samples to collect. }
  \item{...}{currently no additional arguments.  }
}

\value{
  A \pkg{CODA} \code{mcmc} matrix object with columns corresponding to each
  parameter and posterior samples held in the rows.
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

##Now with bayes.lm.ref
n.samples <- 500

out.bayes.lm.ref <- bayes.lm.ref(lm.obj, n.samples)

summary(out.bayes.lm.ref)

}

\keyword{misc}
