\name{ggt.sp}
\alias{ggt.sp}
\title{General Gaussian Template for common spatial covariance functions}

\description{
  The function \code{ggt.sp} performs Bayesian analysis of
  Gaussian models with spatially dependent error structure.
}

\usage{
ggt.sp(formula, data = parent.frame(), coords, run.control, beta.control, cov.model, var.update.control,
                     DIC=TRUE, DIC.start = 1, verbose=TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. The
    details of model specification are given below.  }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{ggt.sp} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{run.control}{a list with tag \code{n.samples} the value of which
    is the number of MCMC iterations.  }
  \item{beta.control}{a list with tags
    \code{beta.update.method}, \code{tuning.matrix},
    \code{beta.starting}, \code{beta.prior}, and \code{sample.order}. The value of
    \code{beta.update.method} is either \code{"gibbs"} or
    \code{"mh"}. With the \code{"mh"} option, the value of the optional tag \code{tuning.matrix} is
    the upper-triangle of the Cholesky decomposed \eqn{p \times p}{p x p}
    tuning matrix, where \eqn{p} is the number of regression parameters.  The
    tuning matrix controls the proposal step size in the
    Metropolis-Hastings sampling.  The tuning values along the diagonal of this matrix
    correspond to the order of the regressors in the model
    formula. Off diagonal elements are the covariance among the terms.  The optional tag \code{beta.starting}
    receives a real vector of length \eqn{p}, whose elements define the starting values of the regressors used in the MCMC
    sampling. The order of the starting values in this vector corresponds
    to the order of the regressors in the model formula.  The optional
    tag \code{beta.prior} is a valid \code{\link{prior}} object that
    defines a normal distribution.  If \code{beta.prior} is not included
    in the \code{beta.control} list then the prior distribution on the
    regressors is assumed flat (i.e., proportional to 1.0).  With the \code{"mh"} option, the value of the optional tag
    \code{sample.order} determines the order in which parameters are updated by sequential
    Metropolis-Hastings steps, see below for details. }
  \item{cov.model}{a quoted key word that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }  
  \item{var.update.control}{a list with each tag corresponding to a
    parameter's name in the specified \code{cov.model}.  The value
    portion of each of these tags is another list with tags
    \code{sample.order}, \code{starting}, \code{tuning}, and
    \code{prior}.  The value of the optional tag \code{sample.order}
    determines the order in which parameters are updated by sequential
    Metropolis-Hastings steps, see below for details.  The \code{sample.order} value is
    relative to the sample order specified for the other
    parameters in the \code{var.update.control} list. The value of the optional tag
    \code{starting} defines the starting value of the parameter in the
    Metropolis-Hastings sampling.  The value associated with the
    \code{tuning} tag defines the step size of the proposal used in the
    Metropolis-Hastings sampling.  The value of the \code{prior} tag is a
    valid \code{\link{prior}} object which describes the prior distribution of
    the parameter.  Any parameter in the specified
    \code{cov.model} can be fixed by using the tag \code{fixed}
    with a value equal to the desired fixed value.  See example below. }
  \item{DIC}{if \code{TRUE}, DIC and associated statistics will be
    calculated to assess model fit. }
  \item{DIC.start}{the sampler iteration to start the DIC calculation.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{ggt.sp}, which is a list with the following tags:
  \item{cov.model}{the covariance model specified by \code{cov.model}. }
  \item{no.tausq}{a logical value indicating if the parameter
    \eqn{\tau^2}{tau^2} was used.  }
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{X}{the \eqn{n \times p}{n x p} matrix of regressors specified by \code{formula}. }
  \item(Y){the the \eqn{n \times 1}{n x 1} response variable vector specified
    by \code{formula}. }
  \item{p.samples}{a matrix of posterior sample for the regressors and
    parameters used in the specified \code{cov.model}.  }
  \item{acceptance}{a matrix of the Metropolis-Hastings sampling
    acceptance rate.  Row names correspond to the parameters used in the
    \code{var.function} and if the \code{beta.update.method} is set
    to \code{"mh"} the acceptance rate for Metropolis-Hastings block
    proposal of the regressors is included.  }
  \item{DIC}{a matrix that holds \eqn{DIC} and associated values \eqn{\bar{D}}{bar(D)},
    \eqn{D[\bar{\theta}]}{D[bar(theta)]}, and \eqn{pD} (Banerjee et al. 2004
    and Spiegelhalter et al. 2002). }
}

\details{
  The Gaussian model assumes that \eqn{y \sim N(X \beta,
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
    I)}{MVN(0, tau^2 I)}, respectively.  Once the marginalized model is fit with \code{ggt.sp} the
  vector of random spatial effects can be recovered with
  \code{\link{sp.effect}} function (see Banerjee et al. 2004, Chapter 5 for
  further details).

  The parameters associated with the \eqn{X} are regressors and defined in the \code{beta.control} list. The \code{cov.model} defines
  the structure of \eqn{n \times n}{n x n} \eqn{\Sigma_y}.  Covariance models
  and associated parameters are: 
  
  \itemize{
    \item{Exponential: }{\eqn{\sigma^2 exp(-\phi d)+\tau^2}{sigma^2 exp(-phi d)+tau^2}, for \eqn{d > 0}{d > 0},
      where \eqn{sigma^2}{sigma^2} is the
      partial sill, \eqn{\tau^2}{tau^2} is the nugget,
      \eqn{\phi}{phi} is the decay parameter, and \eqn{d} is the distance
      between any two sites. }
    
    \item{Gaussian: }{\eqn{\sigma^2 exp(-\phi^2 d^2)+\tau^2}{sigma^2 exp(-phi^2 d^2) + tau^2}, for \eqn{d > 0}{d > 0},
      where \eqn{sigma^2}{sigma^2} is the
      partial sill, \eqn{\tau^2}{tau^2} is the nugget,
      \eqn{\phi}{phi} is the decay parameter, and \eqn{d} is the distance
      between any two sites. }
    
    \item{Spherical: }{\eqn{\sigma^2 (1 - 3/2 \phi d + 1/2(\phi
	d)^3)+\tau^2}{sigma^2(1 - 3/2 phi d + 1/2(phi d)^3)+tau^2},
      for \eqn{0 < d <= 1/\phi}{0 < d <= 1/phi}, where \eqn{sigma^2}{sigma^2} is the
      partial sill, \eqn{\tau^2}{tau^2} is the nugget,
      \eqn{\phi}{phi} is the decay parameter, and \eqn{d} is the distance
      between any two sites. }
    
    \item{Matern: }{\eqn{\sigma^2/(2^(\nu - 1) \Gamma(\nu)) (2
	\nu^0.5 d \phi)^\nu K_\nu (2 \nu^0.5 d \phi)}{sigma^2/(2^(nu -
	1) gamma(nu)) (2 nu^0.5 d phi)^nu K_nu(2 nu^0.5 d phi)}, for \eqn{d > 0}{d > 0}, where \eqn{sigma^2}{sigma^2} is the
      partial sill, \eqn{\tau^2}{tau^2} is the nugget, \eqn{\phi}{phi}
      is the decay parameter, \eqn{\nu}{nu} is decay smoothness,
      \eqn{\Gamma}{gamma} is the gamma function, and \eqn{K_\nu}{K_nu}
      is the modified Bessel of the second kind. }
  }
  
  For all covariance models the \eqn{\tau^2}{tau^2} parameter is
  optional. To omit \eqn{\tau^2}{tau^2} from the model, do not include it in the \code{var.update.control}
  
  The update method for the regressors is either the
  Metropolis-Hastings algorithm or Gibbs.  If the Metropolis-Hastings
  algorithm is chosen, then updates come from a single block draw from a
  multi-variate normal distribution with step size determined by the
  tuning matrix.  Otherwise, a Gibbs sampler is used to sample from the conditional
  posterior distribution of the regressors given the error terms.

  The Metropolis-Hastings algorithm is used to estimate the joint posterior
  distribution of the model's variance and, if specified, regression parameters.  In the simplest case, for each iteration of the algorithm
  a single a single vector of candidate values is drawn from a multivariate normal
  distribution with variance matrix specified by the tuning values
  in the \code{var.update.control} and \code{beta.control} lists.  This
  single block draw is specified by setting the \code{sample.order}
  values equal in the \code{var.update.control} and \code{beta.control}
  lists.  In so doing, a single acceptance rate is monitored for all
  parameters. An acceptance rate of about 30 percent is generally
  recommended (see Gelman et al., 1996). Hand tuning can be tricky, and it is often useful to update parameters independently or
  in smaller sets (i.e., blocks).  This is accomplished by setting the
  \code{sample.order} accordingly. For instance, to update variance
  parameters as a single block separate from the regressors, set
  the all \code{sample.order} values in the \code{var.update.control}
  list equal and different than the \code{sample.order} in the
  \code{beta.control} list.  In this way, two separate acceptance rates
  are monitored and reported.  The extreme is setting all
  \code{sample.order} values different, which allows each parameter to
  be monitored independently. However, the likelihood is evaluated for
  each each parameter block, this requires an inversion of the
  covariance matrix, which if \eqn{n} is large could take a significant amount of
  time.   
}

\seealso{
  \code{\link{prior}}, \code{\link{sp.effect}},
  \code{\link{sp.predict}}, \code{\link{sp.predict}},
  \code{\link{sp.DIC}}, \code{\link{ggt}}
}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

  Gelman, A., Roberts, G.O., and Gilks, W.R. (1996). Efficient Metropolis
  jumping rules. In Bayesian Statistics 5, eds. J.M. Bernardo,
  J.O. Berger, A.P. Dawid, and A.F.M. Smith. Oxford: Oxfod University
  Press, pp. 599-607.
  
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
##Subset data
###############################################
set.seed(1234)

data(BEF)

n.subset <- 100
BEF.subset <- BEF[sample(1:nrow(BEF),n.subset),]

##############################################
##General ggt.sp setup and model fit
##############################################

##Specify the priors, hyperparameters, and variance parameter starting values.
sigmasq.prior <- prior(dist="IG", shape=2, scale=30)
tausq.prior <- prior(dist="IG", shape=2, scale=30)
##The prior on phi corresponds to a prior of 500-2000 meters
##for the effective range (i.e., -log(0.05)/0.0015, -log(0.05)/0.006, when
##using the exponential covariance function).
phi.prior <- prior(dist="LOGUNIF", a=0.0015, b=0.006)

var.update.control <- list("sigmasq"=list(sample.order=0, starting=30, tuning=0.5, prior=sigmasq.prior),
                           "tausq"=list(sample.order=0, starting=30, tuning=0.5, prior=tausq.prior),
                           "phi"=list(sample.order=0, starting=0.006, tuning=0.8, prior=phi.prior))

##Specify the number of samples.
run.control <- list("n.samples" = 1500)

##Specify the model, assign the prior, and starting values for the regressors.
model <- BE_BA_AREA~ELEV+SPR_02_TC1+SPR_02_TC2+SPT_02_TC3

##Note, precision of 0 is same as flat prior alternatively if you want a flat prior
##do not include the beta.prior tag in the beta.control list.
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

##Get the effective range of spatial dependence.
eff.range <- -log(0.05)/ggt.sp.out$p.samples[,"phi"]

##Plot the chain with coda.
mcmc.obj <- mcmc(cbind(ggt.sp.out$p.samples, eff.range))
plot(mcmc.obj)

##############################################
##ggt.sp setup and model fit with a
##fixed variance parameter
##############################################

var.update.control <- list("sigmasq"=list(sample.order=0, starting=30, tuning=0.5, prior=sigmasq.prior),
                           "tausq"=list(sample.order=0, starting=30, tuning=0.5, prior=tausq.prior),
                           "phi"=list(fixed=0.001))

ggt.sp.out <-ggt.sp(formula=model, data=BEF.subset,
                    coords=model.coords,
                    var.update.control=var.update.control,
                    beta.control=beta.control,
                    cov.model="exponential",
                    run.control=run.control, DIC=TRUE, verbose=TRUE)

print(ggt.sp.out$accept)
print(ggt.sp.out$DIC)

##Summarize and plot chain with coda
mcmc.obj <- mcmc(ggt.sp.out$p.samples)
summary(mcmc.obj)
##plot(mcmc.obj)



##############################################
##ggt.sp setup and model fit with single
##block update of the variance and regression
##parameters
##############################################

var.update.control <- list("sigmasq"=list(sample.order=0, starting=30, tuning=0.5, prior=sigmasq.prior),
                           "tausq"=list(sample.order=0, starting=30, tuning=0.5, prior=tausq.prior),
                           "phi"=list(sample.order=0, starting=0.006, tuning=0.8, prior=phi.prior))

#Using Metropolis-Hastings and default flat prior for regression parameters.
beta.control <- list(beta.update.method="mh", sample.order=0)

ggt.sp.out <-ggt.sp(formula=model, data=BEF.subset,
                    coords=model.coords,
                    var.update.control=var.update.control,
                    beta.control=beta.control,
                    cov.model="exponential",
                    run.control=run.control, DIC=TRUE, verbose=TRUE)

##Note single acceptance rate.
print(ggt.sp.out$accept)

}

\keyword{misc}
