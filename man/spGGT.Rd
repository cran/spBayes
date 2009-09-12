\name{spGGT}
\alias{spGGT}
\title{Function for fitting univariate and multivariate
  Bayesian spatial regression models}

\description{
  The function \code{spGGT} fits Gaussian univariate and multivariate Bayesian spatial regression models.
}

\usage{
spGGT(formula, data = parent.frame(), coords, run.control,
       var.update.control, beta.update.control, cov.model,
       verbose=TRUE, ...)
}

\arguments{
  \item{formula}{for a univariate model, this is a symbolic description of the regression model to be
    fit. For a multivariate model, this is a list of univariate formulas.
    See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spGGT} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{run.control}{a list with tags: \code{n.samples} the value of which
    is the number of MCMC iterations; \code{sp.effects} a logical value
    indicating if spatial random effects should be recovered; \code{DIC}
    a logical value indicating if marginalized DIC and associated
    statistics should be computed; \code{DIC.start} and \code{DIC.end}
    the MCMC iteration interval over which DIC should be calculated; and \code{linpack} a logical value
    indicating if LINPACK or LAPACK should be used for matrix algebra
    operations.  The default is to use LINPACK. }
  \item{var.update.control}{a list with each tag corresponding to a
    parameter name used to define the variance structure of the
    model. Valid list tags are \code{K}, \code{Psi}, \code{phi}, and
    \code{nu}.  The value portion of each of these tags is another list with tags
    \code{sample.order}, \code{starting}, \code{tuning}, and
    \code{prior}.  The positive scalar value of the optional tag \code{sample.order}
    determines the order in which parameters are updated by sequential
    Metropolis-Hastings steps.  The \code{sample.order} value is
    relative to the sample order specified for the other
    parameters. The value of the \code{starting} tag defines the starting value of the parameter in the
    Metropolis-Hastings sampling.  If the list element defines a vector or matrix
    of parameters, then the starting value should be set
    accordingly.  If a scalar starting value is given for a vector or
    matrix of parameters then the scalar value is recycled.  If the
    model defines \code{K} or \code{Psi} as a full cross-covariance
    matrix, opposed to a diagonal matrix or scalar, then the starting value should
    be the lower-triangle of the Cholesky square root of the desired starting
    matrix value.  If  \code{K} or \code{Psi} is a diagonal matrix,
    then the starting values should be passed as a vector or scalar.  The value associated with the
    \code{tuning} tag defines the step size of the proposal used in the
    Metropolis-Hastings sampling.  Again the tuning value is recycled if
    needed.  If the parameter element is a matrix, i.e., \code{K} or \code{Psi}, then
    the diagonal elements in the tuning matrix correspond to the column
    major lower-triangle elements in the parameter matrix.  The value of the \code{prior} tag is a
    valid \code{\link{prior}} object which describes the prior distribution of
    the parameter.  Any parameter in the specified model can
    be fixed by using the \code{FIXED} \code{\link{prior}} with the \code{starting} value set to the desired fixed value.  See example below. }
  \item{beta.update.control}{a list with tags
    \code{update}, \code{tuning},
    \code{starting}, \code{prior}, and \code{sample.order}. The value of
    \code{update} is either \code{"GIBBS"} or
    \code{"MH"}. With the \code{"MH"} option, the value of the optional tag \code{tuning} is
    the upper-triangle of the Cholesky decomposed \eqn{p \times p}{p x p}
    tuning matrix, where \eqn{p} is the number of regression parameters.  The
    tuning matrix controls the proposal step size in the
    Metropolis-Hastings sampling.  The tuning values along the diagonal of this matrix
    correspond to the order of the regressors in the model
    formula. Off-diagonal elements are the covariance among the terms.  The optional tag \code{starting}
    receives a vector of length \eqn{p}, whose elements define the starting values of the regressors used in the MCMC
    sampling. The order of the starting values in this vector corresponds
    to the order of the regressors in the model formula.  The optional
    tag \code{prior} is a valid \code{\link{prior}} object that
    defines either a normal or flat distribution.  If \code{prior} is not included
    in the \code{beta.update.control} list then the prior distribution on the
    regressors is assumed flat (i.e., proportional to 1.0).  With the \code{"MH"} option, the value of the optional tag
    \code{sample.order} determines the order in which parameters are updated by sequential
    Metropolis-Hastings steps, see below for details. }
  \item{cov.model}{a quoted key word that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }  
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{spGGT}, which is a list with the following tags:
  \item{cov.model}{the covariance model specified by \code{cov.model}. }
  \item{no.Psi}{a logical value indicating if the parameter
    \code{Psi} was used.  }
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{X}{the \eqn{n \times p}{n x p} matrix of regressors specified by \code{formula}. }
  \item{Y}{the the \eqn{n \times 1}{n x 1} response variable vector specified
    by \code{formula}. }
  \item{p.samples}{a matrix of posterior samples for the defined
    parameters.  If \code{K} or \code{Psi} is specified as a full
    cross-covariance matrix then the corresponding samples in
    \code{p.samples} are column major lower-triangle elements of the
    matrix. }
  \item{acceptance}{a matrix of the Metropolis-Hastings sampling
    acceptance rate.  Row names correspond to the parameters updated with Metropolis-Hastings. }
  \item{DIC}{a matrix that holds marginalized DIC and associated values
    is returned if \code{DIC} is true in the \code{run.control} list. }
  \item{sp.effects}{a matrix that holds samples from the posterior
    distribution of the response specific spatial random effects. The rows of this matrix
    correspond to the \eqn{n}{n} point observations and the columns are the
    posterior samples.  If a multivariate model is fit to \eqn{m}{m}
    response variables, the \code{sp.effects} matrix has
    \eqn{mn}{mn} rows. The random effects samples for points are held in rows
    1:m, (m+1):2m, ..., ((i-1)m+1):im, ..., ((n-1)m+1):nm, where i
    = 1 ... n (e.g., the samples for the first point are in rows 1:m,
    second point in rows (m+1):2m, etc.). }


}

\details{
  Please refer to Section 3.1 in the vignette.
}

\seealso{
  \code{\link{prior}}, \code{\link{spPredict}}, \code{\link{spDIC}}
}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

  %Further information on the package \pkg{spBayes} can be found at: \url{http://blue.fr.umn.edu/spatialBayes}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Bradley P. Carlin \email{brad@biostat.umn.edu}.
}

\examples{
\dontrun{

data(FBC07.dat)
Y.1 <- FBC07.dat[1:100,"Y.1"]
Y.2 <- FBC07.dat[1:100,"Y.2"]
coords <- as.matrix(FBC07.dat[1:100,c("coord.X", "coord.Y")])

#############################
##   Univariate models
#############################

##Fit some competing models.
##Full model with nugget (Psi), partial sill (K),
##and spatial range (Phi).
K.prior <- prior(dist="IG", shape=2, scale=5)
Psi.prior <- prior(dist="IG", shape=2, scale=5)
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.5, prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=1000, "sp.effects"=TRUE)

Fit.m1 <- spGGT(formula=Y.2~1, run.control=run.control,
                 coords=coords, var.update.control=var.update.control,
                 beta.update.control=beta.control,
                 cov.model="exponential")

plot(Fit.m1$p.samples)
summary(Fit.m1$p.samples)

##First reduced model with no nugget just partial sill and range.
var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

Fit.m2 <- spGGT(formula=Y.2~1, run.control=run.control,
                 coords=coords, var.update.control=var.update.control,
                 beta.update.control=beta.control,
                 cov.model="exponential")

plot(Fit.m2$p.samples)

##Second reduced model with no nugget just partial sill and fixed
##range at the WRONG value, which forces the wrong estimate of K.
phi.prior <- prior(dist="FIXED")
       
var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "phi"=list(starting=0.1, prior=phi.prior)
       )

Fit.m3 <- spGGT(formula=Y.2~1, run.control=run.control,
                 coords=coords, var.update.control=var.update.control,
                 beta.update.control=beta.control,
                 cov.model="exponential")
plot(Fit.m3$p.samples)

##For this data, DIC cannot distinguish between the first two
##models but does show that the third model is clearly
##wrong.  An empirical semivariogram is useful
##for recognizing that the full model is the correct choice
##over the second model.
print(spDIC(Fit.m1, DIC.unmarg=FALSE))
print(spDIC(Fit.m2, DIC.unmarg=FALSE))
print(spDIC(Fit.m3, DIC.unmarg=FALSE))

##Use akima or MBA to produce interpolated surfaces of
##estimated random spatial effects.

m1.surf <- mba.surf(cbind(coords, rowMeans(Fit.m1$sp.effects)),
                    no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(m1.surf, xaxs="r", yaxs="r",
      main="Y.2 random spatial effects")
points(coords, pch=19)

##Now the full model using the Matern spatial correlation function.
##Also, for fun, update the spatial range parameters
##phi and nu in a separate MH block.  This provides finer control
##on the parameters' acceptance rate, but takes a bit longer to
##run.
phi.prior <- prior(dist="UNIF", a=0.06, b=3)
nu.prior <- prior(dist="UNIF", a=0.01, b=2)

var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.5, prior=Psi.prior),
       "phi"=list(sample.order=1, starting=0.1,
         tuning=0.5, prior=phi.prior),
       "nu"=list(sample.order=1, starting=0.1,
         tuning=0.5, prior=nu.prior)     
       )

run.control <- list("n.samples"=1000, "sp.effects"=FALSE)

Fit.m4 <-
  spGGT(formula=Y.2~1, run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="matern")

plot(Fit.m4$p.samples)

##Solve the Matern correlation function for the estimated
##95\% CI for effective range (e.g., distance at which the
##correlation drops to 0.05).
phi.hat <-  quantile(Fit.m4$p.samples[,"Phi"], c(0.50,0.025,0.975))
nu.hat <- quantile(Fit.m4$p.samples[,"Nu"], c(0.50,0.025,0.975))

matern <- function(cor, eff.d, phi, nu){
  cor - (eff.d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=eff.d*phi, nu=nu)
}

##50% (2.5%, 97.5%)
print(uniroot(matern, interval=c(1, 50),
              cor=0.05, phi=phi.hat[1], nu=nu.hat[1])$root)
print(uniroot(matern, interval=c(1, 50),
              cor=0.05, phi=phi.hat[2], nu=nu.hat[2])$root)
print(uniroot(matern, interval=c(1, 50),
              cor=0.05, phi=phi.hat[3], nu=nu.hat[3])$root)

##Finally, the full model again but this time update the Beta
##using MH.  Using the default of sample.order=0, this sampling
##scheme is simply a single block MH proposal of all model parameters.
##Also, the first run uses the default tuning value for the Beta
##(i.e., chol(vcov(lm(Y~X))). Both chains use a normal prior on Beta;
##however, the first chain allows for a vague prior, where as the
##second is much too strict.
var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.5, prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

##Chain 1: vague prior normal variance of 1/0.001 = 1000
beta.prior <- prior(dist="NORMAL", mu=5, precision=0.001)
beta.control <- list(update="MH", prior=beta.prior)

Fit.m5.a <-
  spGGT(formula=Y.2~1, run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

##Chain 2: restrictive prior normal variance of 1/5 = 0.2
##(i.e., crazy strict)
beta.prior <- prior(dist="NORMAL", mu=5, precision=5)
beta.control <- list(update="MH", prior=beta.prior, tuning=0.75)

Fit.m5.b <-
  spGGT(formula=Y.2~1, run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

plot(mcmc.list(Fit.m5.a$p.samples, Fit.m5.b$p.samples))


##Now just for fun, here is the full model but using the weakly
##informative half-Cauchy for the nugget (Psi) and partial sill (K).

K.prior <- prior(dist="HC", a=50)
Psi.prior <- prior(dist="HC", a=50)
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, sample.order=1, tuning=2, prior=Psi.prior),
       "phi"=list(starting=0.1, sample.order=2, tuning=1, prior=phi.prior)
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=3000, "sp.effects"=TRUE)

Fit.m1 <- spGGT(formula=Y.2~1, run.control=run.control,
                 coords=coords, var.update.control=var.update.control,
                 beta.update.control=beta.control,
                 cov.model="exponential")

plot(Fit.m1$p.samples)
summary(Fit.m1$p.samples)


#############################
##   Multivariate models
#############################

##Full model with non-spatial cross-covariance
##(Psi) matrix, spatial cross-covariance matrix (K),
##and response specific spatial range parameter (Phi)
##(i.e., a non-separable model).

##These chains really need to be run for several thousand
##samples (e.g., 10,000). Also note that it is much easier to maintain
##a healthy acceptance rate for each parameter if they
##are updated individually using the sample.order
##directive, as shown in the second example below.
K.prior <- prior(dist="IWISH", df=2, S=diag(c(3, 6)))
Psi.prior <- prior(dist="IWISH", df=2, S=diag(c(7, 5)))
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

##The matrix starting values for K and Psi are actually
##passed as the sqrt of the desired starting values.
##This is only done if K or Psi are full cross-covariance
##matrices and serve to insure that the true starting value
##matrix is positive definite.
K.starting <- matrix(c(2,-3, 0, 1), 2, 2)
Psi.starting <- diag(c(3, 2))

var.update.control <-
  list("K"=list(starting=K.starting, tuning=diag(c(0.08, 0.3, 0.08)),
         prior=K.prior),
       "Psi"=list(starting=Psi.starting, tuning=diag(c(0.08, 0.3, 0.08)),
         prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5,
         prior=list(phi.prior, phi.prior))
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=2000, "sp.effects"=FALSE)

Fit.mv.m1 <-
  spGGT(formula=list(Y.1~1, Y.2~1), run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

plot(Fit.mv.m1$p.samples)

##Calculate the mean K correlation matrix.  Note that since K is
##specified as a full cross-covariance matrix, the posterior
##samples are returned for the lower triangle of the symmetric K,
##in column major order.
K <- matrix(0, 2, 2)
K[lower.tri(K, diag=TRUE)] <-
  colMeans(Fit.mv.m1$p.samples[,c("K_1","K_2","K_3")])
K[upper.tri(K, diag=FALSE)] <- t(K)[upper.tri(K, diag=FALSE)]
print(cov2cor(K))


##Now using sample order.  Note this takes ~3 times as long
##to run, but provides much finer control of acceptance rates.
var.update.control <-
  list("K"=list(sample.order=0, starting=K.starting,
         tuning=diag(c(0.3, 0.5, 0.2)), prior=K.prior),
       "Psi"=list(sample.order=1, starting=Psi.starting,
         tuning=diag(c(0.2, 0.5, 0.2)), prior=Psi.prior),
       "phi"=list(sample.order=2, starting=0.1,
         tuning=0.5, prior=list(phi.prior, phi.prior))
       )

run.control <- list("n.samples"=2000, "sp.effects"=TRUE)

Fit.mv.m2 <-
  spGGT(formula=list(Y.1~1, Y.2~1), run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

plot(Fit.mv.m2$p.samples)

##Now again check out the random spatial effects.  Recall,
##these are stored in the sp.effects matrix, which is organized
##as m consecutive rows for each point.
rand.effect.Y.1 <-
  Fit.mv.m2$sp.effects[seq(1, nrow(Fit.mv.m2$sp.effects), 2),]

rand.effect.Y.2 <-
  Fit.mv.m2$sp.effects[seq(2, nrow(Fit.mv.m2$sp.effects), 2),]

rand.effect.Y.1.surf <-
  mba.surf(cbind(coords, rowMeans(rand.effect.Y.1)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

rand.effect.Y.2.surf <-
  mba.surf(cbind(coords, rowMeans(rand.effect.Y.2)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

##The estimated negative value of off-diagonal element in K is
##clearly seen in the map of random spatial effects.
par(mfrow=c(1,2))
image(rand.effect.Y.1.surf, xaxs="r", yaxs="r",
      main="Y.1 random spatial effects")
contour(rand.effect.Y.1.surf, add=TRUE)
image(rand.effect.Y.2.surf, xaxs="r", yaxs="r",
      main="Y.2 random spatial effects")
contour(rand.effect.Y.2.surf, add=TRUE)

##Now a simpler model.  Specifically, the full model shows that
##there is negligible non-spatial cross-covariance, therefore,
##Psi might just be modeled with a diagonal covariance matrix.
##This is actually the model which gave rise to this synthetic
##data set.
K.prior <- prior(dist="IWISH", df=2, S=diag(c(3, 6)))
Psi.prior <- prior(dist="IG", shape=2, scale=5)
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(sample.order=0, starting=K.starting,
         tuning=diag(c(0.3, 0.5, 0.2)), prior=K.prior),
       "Psi"=list(sample.order=1, starting=5,
         tuning=0.5, prior=list(Psi.prior, Psi.prior)),
       "phi"=list(sample.order=2, starting=0.1,
         tuning=0.5, prior=list(phi.prior, phi.prior))
       )

run.control <- list("n.samples"=2000)

Fit.mv.m3 <-
  spGGT(formula=list(Y.1~1, Y.2~1), run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

plot(Fit.mv.m3$p.samples)
summary(Fit.mv.m3$p.samples)

##Finally, a really simple model which we know makes no sense
##for this data set.  Specifically, a single variance parameter
##for the diagonal elements of K and Psi and common spatial range,
##phi.  Also, we use a spherical spatial covariance function, just
##because we can.
K.prior <- prior(dist="IG", shape=2, scale=5)
Psi.prior <- prior(dist="IG", shape=2, scale=5)
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.5, prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

Fit.mv.m4 <-
  spGGT(formula=list(Y.1~1, Y.2~1), run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="spherical")

plot(Fit.mv.m4$p.samples)

}
}

\keyword{misc}
