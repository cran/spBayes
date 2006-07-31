\name{bayes.lm.ref}
\alias{bayes.lm.ref}          
\title{Simple Bayesian linear model with non-informative priors }

\description{
  Given an \code{lm} object, the \code{bayes.lm.ref} function fits a
  simple Bayesian linear model with reference (non-informative) priors.  
}

\usage{
  bayes.lm.ref(lm.obj, n.samples)
}

\arguments{
  \item{lm.obj}{an object returned by \code{lm}. }
  \item{n.samples}{the number of posterior samples to collect. }
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
  Andrew O. Finley \email{afinley@stat.umn.edu}, \cr
  Bradley P. Carlin \email{brad@biostat.umn.edu}.
}

\examples{
## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
summary(lm.D90 <- lm(weight ~ group - 1))

##Now with bayes.lm.ref
fit <- bayes.lm.ref(lm.D90, 1000)

summary(fit)
plot(fit)
}

\keyword{misc}
