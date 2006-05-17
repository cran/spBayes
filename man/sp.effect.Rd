\name{sp.effect}
\alias{sp.effect}
\title{Recovers the random spatial effects from a ggt.sp object}

\description{
  The function \code{sp.effect} recovers the vector of random
  spatial effects from a given marginalized model contained in \code{\link{ggt.sp}} object.
}

\usage{
sp.effect(ggt.sp.obj, start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{ggt.sp.obj}{an object returned by \code{\link{ggt.sp}} (i.e.,
    object of class \code{ggt.sp}). }
  \item{start}{specifies the first sample included in spatial process recovery.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{end}{specifies the last sample included in spatial process recovery.
    The default is to include from \code{start} to
    \code{nrow(ggt.sp.obj$p.samples)}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10}
    then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{verbose}{if \code{TRUE} calculation progress is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{...}{currently no additional arguments.  }
}

\value{
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix of the observation
    coordinates from \code{ggt.sp.obj$coords}.  }
  \item{pp.samples}{an \eqn{n \times m}{n x m} matrix that holds \eqn{m}
    samples from the posterior predictive distribution of the spatial
    process. See below for details. }
}

\details{
  The Gaussian model fit by \code{\link{ggt.sp}} assumes that \eqn{y \sim N(X \beta,
    \Sigma_y)}{y~N(X \beta, \Sigma_y)},
  where \eqn{y} is the \eqn{n \times 1}{n x 1}
  response vector, \eqn{X} is the \eqn{n \times p}{n x p}
  matrix of regressors, and \eqn{\Sigma_y = (\sigma^2
    R(\phi)+\tau^2 I)}{\Sigma_y=(sigma^2 R(phi)+tau^2 I)}, where \eqn{R} is a
  valid correlation function specified by \code{cov.model}. Under this
  marginalized model, we assume \eqn{\Sigma_y
    \sim MVN(0, \sigma^2 R(\phi)+\tau^2 I)}{\Sigma_y ~ MVN(0,
    sigma^2 R(phi)+tau^2 I)}.  The unmarginalized model partitions the
  variance terms into spatial and non-spatial components which follow
  \eqn{MVN(0, \sigma^2 R(\phi))}{MVN(0, sigma^2 R(phi))} and \eqn{MVN(0, \tau^2
    I)}{MVN(0, tau^2 I)}, respectively.  Once the marginalized model is fit with \code{\link{ggt.sp}} the
  vector of random spatial effects can be recovered with \code{sp.effect} function.
}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
  
  Spiegelhalter, D.J., Best, N., Carlin, B.P., and van der Linde, A. (2002). Bayesian measures of model complexity and fit (with discussion). J. Roy. Statist. Soc., Ser. B, 64, 583-639.

  Further information on the package \pkg{spBayes} can be found
  at: \url{http://blue.fr.umn.edu/spatialBayes}.
}

\author{
  Andrew O. Finley \email{afinley@stat.umn.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Bradley P. Carlin \email{brad@biostat.umn.edu}.
}

\examples{

###############################################
##subset data
###############################################
set.seed(1234)

data(BEF)

n.subset <- 100
BEF.subset <- BEF[sample(1:nrow(BEF),n.subset),]

##############################################
##general ggt.sp setup and model fit
##############################################

##Specify the priors, hyperparameters, and variance parameter starting values.
sigmasq.prior <- prior(dist="IG", shape=2, scale=30)
tausq.prior <- prior(dist="IG", shape=2, scale=30)
##the prior on phi corresponds to a prior of 500-2000 meters
##for the effective range (i.e., -log(0.05)/0.0015, -log(0.05)/0.006, when
##using the exponential covariance function).
phi.prior <- prior(dist="LOGUNIF", a=0.0015, b=0.006)

var.update.control <- list("sigmasq"=list(sample.order=0, starting=30, tuning=0.5, prior=sigmasq.prior),
                           "tausq"=list(sample.order=0, starting=30, tuning=0.5, prior=tausq.prior),
                           "phi"=list(sample.order=0, starting=0.006, tuning=0.8, prior=phi.prior))

##specify the number of samples
run.control <- list("n.samples" = 1500)

##specify some model, assign the prior and starting values for the regressors
model <- BE_BA_AREA~ELEV+SPR_02_TC1+SPR_02_TC2+SPT_02_TC3

##note, precision of 0 is same as flat prior.
beta.prior <- prior(dist="NORMAL", mu=rep(0, 5), precision=diag(0, 5))

beta.control <- list(beta.update.method="gibbs", beta.starting=rep(0, 5), beta.prior=beta.prior)

model.coords <- cbind(BEF.subset$XUTM, BEF.subset$YUTM)

ggt.sp.out <-ggt.sp(formula=model, data=BEF.subset,
                    coords=model.coords,
                    var.update.control=var.update.control,
                    beta.control=beta.control,
                    cov.model="exponential",
                    run.control=run.control, DIC=TRUE, verbose=TRUE)

print(ggt.sp.out$accept)
print(ggt.sp.out$DIC)

##############################################
##recover random spatial effect
##Y = X^tB + w + e, where w is the random
##spatial effect and e is pure error
##############################################

##get the spatial random effect w
sp.effect.out <- sp.effect(ggt.sp.out, start=500, thin=2)
w <- rowMeans(sp.effect.out$pp.samples)

par(mfrow=c(1,2))
int.obs <- interp(BEF.subset$XUTM, BEF.subset$YUTM, BEF.subset$BE_BA_AREA)
image(int.obs, xlab="Longitude (meters)", ylab="Latitude (meters)", main="Observed density\nof American beech (Y)")
contour(int.obs, add=TRUE)

int.w <- interp(BEF.subset$XUTM, BEF.subset$YUTM, w)
image(int.w, xlab="Longitude (meters)", ylab="Latitude (meters)", main="Spatial random effect density (w)")
contour(int.w, add=TRUE)

}

\keyword{misc}
