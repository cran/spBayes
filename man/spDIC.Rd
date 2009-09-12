\name{spDIC}
\alias{spDIC}          
\title{Calculates DIC and associated statistics given a model object}

\description{
  The function \code{spDIC} calculates model DIC and associated
  statistics given a \code{\link{bayesLMRef}}, \code{\link{bayesLMConjugate}}, \code{\link{spGGT}}, \code{\link{spLM}},  
  \code{\link{spMvLM}}, \code{\link{bayesGeostatExact}}, or
  \code{\link{spGLM}} object.
}

\usage{
  spDIC(sp.obj, start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{bayesLMRef}}, \code{\link{bayesLMConjugate}}, \code{\link{spGGT}},
    \code{\link{spLM}}, \code{\link{spMvLM}},
  \code{\link{bayesGeostatExact}}, or \code{\link{spGLM}}}
  \item{start}{specifies the first sample included in the DIC calculation.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{end}{specifies the last sample included  in the prediction calculation.
    The default is to use all posterior samples in \code{sp.obj}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10}
    then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{verbose}{if \code{TRUE} calculation progress is printed to the
    screen; otherwise, nothing is printed to the screen.  }
  \item{...}{currently no additional arguments.  }
}

\value{
  A list with some of the following tags:
  \item{DIC}{a matrix holding DIC and associated statistics}
  \item{sp.effects}{if \code{sp.obj}
  specifies a spatial model without pre-calculated spatial effects then \code{spDIC} calculates the
    spatial effects. }
}

\details{
  Please refer to Section 3.3 in the vignette.
}

\seealso{\code{\link{bayesLMRef}}, \code{\link{bayesLMConjugate}}, \code{\link{spGGT}}, \code{\link{bayesGeostatExact}}, \code{\link{spLM}}, \code{\link{spGLM}}}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}.
}

\examples{
\dontrun{
###########################################
##          DIC for spLM
###########################################

##Use some more observations
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])

##############################
##With and without predictive
##process
##############################
m.1a <- spLM(Y~1, coords=coords, knots=c(5, 5, 0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=FALSE, n.report=100)

print(spDIC(m.1a, start=100, thin=2)$DIC)

m.1b <- spLM(Y~1, coords=coords, knots=c(7, 7, 0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=FALSE, n.report=100)

print(spDIC(m.1b, start=100, thin=2)$DIC)

m.2 <- spLM(Y~1, coords=coords, 
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=FALSE, n.report=100)

print(spDIC(m.2, start=100, thin=2)$DIC)

###########################################
##    DIC for bayesGeostatExact
###########################################
data(FORMGMT.dat)

n = nrow(FORMGMT.dat)
p = 5 ##an intercept an four covariates

n.samples <- 10

coords <- cbind(FORMGMT.dat$Longi, FORMGMT.dat$Lat)

phi <- 3/0.07

beta.prior.mean <- rep(0, times=p)
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

alpha <- 1/1.5

sigma.sq.prior.shape <- 2.0
sigma.sq.prior.rate <- 10.0

##With covariates
m.3 <-
  bayesGeostatExact(Y~X1+X2+X3+X4, data=FORMGMT.dat,
                      n.samples=n.samples,
                      beta.prior.mean=beta.prior.mean,
                      beta.prior.precision=beta.prior.precision,
                      coords=coords, phi=phi, alpha=alpha,
                      sigma.sq.prior.shape=sigma.sq.prior.shape,
                      sigma.sq.prior.rate=sigma.sq.prior.rate,
                      sp.effects=FALSE, verbos=FALSE)

print(spDIC(m.3)$DIC)

##Without covariates
p <- 1 ##intercept only
beta.prior.mean <- 0
beta.prior.precision <- 0

m.4 <-
  bayesGeostatExact(Y~1, data=FORMGMT.dat,
                      n.samples=n.samples,
                      beta.prior.mean=beta.prior.mean,
                      beta.prior.precision=beta.prior.precision,
                      coords=coords, phi=phi, alpha=alpha,
                      sigma.sq.prior.shape=sigma.sq.prior.shape,
                      sigma.sq.prior.rate=sigma.sq.prior.rate,
                      sp.effects=FALSE, verbose=FALSE)

print(spDIC(m.4)$DIC)

##Lower DIC is better, so go with the covariates.

###########################################
##         DIC for spGGT
###########################################
data(FBC07.dat)

Y.2 <- FBC07.dat[1:100,"Y.2"]
coords <- as.matrix(FBC07.dat[1:100,c("coord.X", "coord.Y")])

##m.5 some model with spGGT.
K.prior <- prior(dist="IG", shape=2, scale=5)
Psi.prior <- prior(dist="IG", shape=2, scale=5)
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(starting=5, tuning=0.1, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.1, prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=1000, "sp.effects"=TRUE)

m.5 <-
  spGGT(formula=Y.2~1, run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential", verbose=FALSE)

DIC <- spDIC(m.5$DIC)
print(DIC)

###########################################
##     Compare DIC between non-spatial
##           and spatial models
###########################################
data(FBC07.dat)
Y.2 <- FBC07.dat[1:150,"Y.2"]
coords <- as.matrix(FBC07.dat[1:150,c("coord.X", "coord.Y")])

##############################
##Non-spatial model
##############################
m.1 <- bayesLMConjugate(Y.2~1, n.samples = 2000,
                          beta.prior.mean=0, beta.prior.precision=0,
                          prior.shape=-0.5, prior.rate=0, verbose=FALSE)

summary(m.1$p.samples)

dic.m1 <- spDIC(m.1)
print(dic.m1)

##############################
##Spatial model
##############################
m.2 <- spLM(Y.2~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.1,"sigma.sq"=5, "tau.sq"=5),
             sp.tuning=list("phi"=0.03, "sigma.sq"=0.03, "tau.sq"=0.03),
             priors=list("phi.Unif"=c(0.06, 3), "sigma.sq.IG"=c(2, 5),
               "tau.sq.IG"=c(2, 5)),
             cov.model="exponential",
             n.samples=2000, verbose=FALSE, n.report=100)

summary(m.2$p.samples)

dic.m2 <- spDIC(m.2)
print(dic.m2)
}
}

\keyword{misc}
