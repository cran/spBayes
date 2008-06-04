\name{sp.DIC}
\alias{sp.DIC}          
\title{Calculates model DIC and associated statistics given a bayes.lm.ref, bayes.lm.conjugate, ggt.sp,
  sp.lm, or bayes.geostat.exact object}

\description{
  The function \code{sp.DIC} calculates model DIC and associated
  statistics given a \code{\link{bayes.lm.ref}}, \code{\link{bayes.lm.conjugate}}, \code{\link{ggt.sp}}, \code{\link{sp.lm}}, or 
  \code{\link{bayes.geostat.exact}} object.
}

\usage{
  sp.DIC(sp.obj, DIC.marg=TRUE, DIC.unmarg=TRUE,
         start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{bayes.lm.ref}}, \code{\link{bayes.lm.conjugate}}, \code{\link{ggt.sp}},
  \code{\link{bayes.geostat.exact}}, or \code{\link{sp.lm}}}
  \item{DIC.marg}{a logical value indicating if marginalized DIC and
    associated statistics should be calculated. Note, this argument is ignored
    when \code{sp.obj} specifies a non-spatial model. }
  \item{DIC.unmarg}{a logical value indicating if unmarginalized DIC and
    associated statistics should be calculated. Note, this argument is ignored
    when \code{sp.obj} specifies a non-spatial model.}
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
  \item{DIC}{a matrix holding DIC and associated statistics when
  \code{sp.obj} specifies a non-spatial model.  }
  \item{DIC.marg}{a matrix holding marginalized DIC and associated
    statistics.  }
  \item{DIC.unmarg}{a matrix holding unmarginalized DIC and associated
    statistics.  }
  \item{sp.effects}{if DIC.ummarg is true and if \code{sp.obj}
  specifies a spatial model without pre-calculated spatial effects then \code{sp.DIC} calculates the
    spatial effects. }
}

\details{
  Please refer to Section 3.3 in the vignette.
}

\seealso{\code{\link{bayes.lm.ref}}, \code{\link{bayes.lm.conjugate}}, \code{\link{ggt.sp}}, \code{\link{bayes.geostat.exact}}, \code{\link{sp.lm}}}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
  
  %Further information on the package \pkg{spBayes} can be found
  %at: \url{http://blue.fr.umn.edu/spatialBayes}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}.
}

\examples{
\dontrun{
###########################################
##          DIC for sp.lm
###########################################

#############################
##Modified predictive process
##############################
##Use some more observations
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])

##############################
##Using unmarginalized DIC
##to assess number of knots
##############################
m.1 <- sp.lm(Y~1, coords=coords, knots=c(5, 5, 0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100, sp.effects=TRUE)

print(sp.DIC(m.1, start=100, thin=2, DIC.marg=TRUE, DIC.unmarg=TRUE))

m.2 <- sp.lm(Y~1, coords=coords, knots=c(7, 7, 0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100, sp.effects=TRUE)

print(sp.DIC(m.2, start=100, thin=2, DIC.marg=TRUE, DIC.unmarg=TRUE))

###########################################
##    DIC for bayes.geostat.exact
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
  bayes.geostat.exact(Y~X1+X2+X3+X4, data=FORMGMT.dat,
                      n.samples=n.samples,
                      beta.prior.mean=beta.prior.mean,
                      beta.prior.precision=beta.prior.precision,
                      coords=coords, phi=phi, alpha=alpha,
                      sigma.sq.prior.shape=sigma.sq.prior.shape,
                      sigma.sq.prior.rate=sigma.sq.prior.rate,
                      sp.effects=FALSE)

print(sp.DIC(m.3, DIC.marg=TRUE, DIC.unmarg=FALSE))

##Without covariates
p <- 1 ##intercept only
beta.prior.mean <- 0
beta.prior.precision <- 0

m.4 <-
  bayes.geostat.exact(Y~1, data=FORMGMT.dat,
                      n.samples=n.samples,
                      beta.prior.mean=beta.prior.mean,
                      beta.prior.precision=beta.prior.precision,
                      coords=coords, phi=phi, alpha=alpha,
                      sigma.sq.prior.shape=sigma.sq.prior.shape,
                      sigma.sq.prior.rate=sigma.sq.prior.rate,
                      sp.effects=FALSE)

print(sp.DIC(m.4, DIC.marg=TRUE, DIC.unmarg=FALSE))

##Lower DIC is better, so go with the covariates.

###########################################
##         DIC for ggt.sp
###########################################
data(FBC07.dat)

Y.2 <- FBC07.dat[1:100,"Y.2"]
coords <- as.matrix(FBC07.dat[1:100,c("coord.X", "coord.Y")])

##m.5 some model with ggt.sp.
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
  ggt.sp(formula=Y.2~1, run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

##Now with the ggt.sp object, m.5, calculate the DIC
##for both the unmarginalized and marginalized models.
##The likelihoods for these models are given by equation 6 and 7
##within the vignette.

DIC <- sp.DIC(m.5)
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
m.1 <- bayes.lm.conjugate(Y.2~1, n.samples = 2000,
                          beta.prior.mean=0, beta.prior.precision=0,
                          prior.shape=-0.5, prior.rate=0)

summary(m.1$p.samples)

dic.m1 <- sp.DIC(m.1)

##> dic.m1
##                  [,1]
##bar.D       503.023678
##D.bar.Omega 501.059965
##pD            1.963713
##DIC         504.987392


##############################
##Spatial model
##############################
m.2 <- sp.lm(Y.2~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.1,"sigma.sq"=5, "tau.sq"=5),
             sp.tuning=list("phi"=0.03, "sigma.sq"=0.03, "tau.sq"=0.03),
             priors=list("phi.Unif"=c(0.06, 3), "sigma.sq.IG"=c(2, 5),
               "tau.sq.IG"=c(2, 5)),
             cov.model="exponential",
             n.samples=2000, verbose=TRUE, n.report=100, sp.effects=TRUE)

summary(m.2$p.samples)

dic.m2 <- sp.DIC(m.2)

##> dic.m2
##$DIC.marg
##                 value
##bar.D       479.904433
##D.bar.Omega 476.856815
##pD            3.047617
##DIC         482.952050

##$DIC.unmarg
##                value
##bar.D       330.50408
##D.bar.Omega 258.64266
##pD           71.86142
##DIC         402.36550


}
}

\keyword{misc}
