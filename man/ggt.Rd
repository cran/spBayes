\name{ggt}
\alias{ggt}
\title{General Gaussian Template}

\description{
  The function \code{ggt} performs Bayesian analysis of Gaussian models with potentially complex hierarchical error structures.
}

\usage{
ggt(formula, data = parent.frame(), run.control, beta.control,
    var.function, var.update.control, DIC=TRUE, DIC.start = 1, verbose=TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression portion of
    the model. The details of model specification are given below.  }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in \code{data}, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{ggt} is called.  }
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
    formula. Off diagonal elements are the covariance among the
    terms.  The optional tag \code{beta.starting}
    receives a real vector of length \eqn{p}, whose elements define the starting values of the regressors used in the MCMC
    sampling. The order of the starting values in this vector corresponds
    to the order of the regressors in the model formula.  The optional
    tag \code{beta.prior} is a valid \code{\link{prior}} object that
    defines a normal distribution.  If \code{beta.prior} is not included
    in the \code{beta.control} list then the prior distribution on the
    regressors is assumed flat (i.e., proportional to 1.0).  With the \code{"mh"} option, the value of the optional tag
    \code{sample.order} determines the order in which parameters are updated by sequential
    Metropolis-Hastings steps, see below for details.}
  \item{var.function}{an R function that defines the error structure
    of the model. The output of this function must be a \eqn{n \times n}{n x n}
    positive definite matrix, where \eqn{n} is the number of
    observations defined by the \code{formula}.  The function must not have any arguments.  Each
    parameter in the function's body must have a corresponding tag in
    the \code{var.update.control} list. Other variables in the function
    body must be available in the environment from which \code{ggt} is called. }
  \item{var.update.control}{a list with each tag corresponding to a
    parameter's name in the \code{var.function}.  The value
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
    the parameter.  }
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
  An object of class \code{ggt}, which is a list with the following tags:
  \item{p.samples}{a matrix of posterior sample for the regressors and
    parameters used in the \code{var.function}.  }
  \item{acceptance}{a matrix of the Metropolis-Hastings sampling
    acceptance rate.  Row names correspond to the parameters used in the
    \code{var.function} and if the \code{beta.update.method} is set
    to \code{"mh"} the acceptance rate for Metropolis-Hastings block
    proposal of the regressors is included.  }
  \item{var.function}{a copy of the \code{var.function}.  }
  \item{DIC}{a matrix that holds \eqn{DIC} and associated values \eqn{\bar{D}}{bar(D)},
    \eqn{D[\bar{\theta}]}{D[bar(theta)]}, and \eqn{pD} (Banerjee et al. 2004
    and Spiegelhalter et al. 2002). }
}

\details{
  The Gaussian model assumes that \eqn{y \sim N(X \beta,
    \Sigma_y)}{y~N(X \beta, \Sigma_y)}, where \eqn{y} is the \eqn{n \times 1}{n x 1}
  response vector, \eqn{X} is the \eqn{n \times p}{n x p} matrix of
  regressors, and \eqn{\Sigma_y = Q_y\sigma^2}. The parameters associated with the
  \eqn{X} are  defined in the \code{beta.control} list.  The \code{var.function} defines
  the structure of \eqn{n \times n}{n x n} \eqn{\Sigma_y}.  For example,
  under the simple linear model, \eqn{Q_y} is the \eqn{n \times n}{n x n}
  identity matrix and \eqn{\sigma^2} is the only parameter defined in
  the \code{var.update.control} list.
  
  The update method for the regressors is either the
  Metropolis-Hastings algorithm or Gibbs.  If the Metropolis-Hastings
  algorithm is chosen, then updates come from a single block draw from a
  multi-variate normal distribution with step size determined by the
  tuning matrix.  Otherwise, a Gibbs sampler is used to sample from the conditional
  posterior distribution of the regressors given the error terms.

  The Metropolis-Hastings algorithm is used to estimate the joint posterior
  distribution of the model's variance and, if specified, regression
  parameters.  In the simplest case, for each iteration of the algorithm
  a single vector of candidate values is drawn from a multivariate normal
  distribution with variance matrix specified by the tuning values
  in the \code{var.update.control} and \code{beta.control} lists.  This
  single block draw is specified by setting the \code{sample.order}
  values equal in the \code{var.update.control} and \code{beta.control}
  lists.  In so doing, a single acceptance rate is monitored for all
  parameters.  An acceptance rate of about 30 percent is generally
  recommended (see Gelman et al., 1996).  Hand tuning can be tricky, and it is often useful to update parameters independently or
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
  \code{\link{ggt.sp}}
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
##General ggt setup and model fit
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

##Get the distance matrix.
D <- as.matrix(dist(cbind(BEF.subset$XUTM, BEF.subset$YUTM)))

##Define the variance function (e.g., exponential covariance)
var.function <- function(){
  sigmasq*exp(-phi*D)+tausq*diag(n.subset)
}

##Specify the number of samples.
run.control <- list("n.samples" = 1500)

##Specify the model, assign the prior, and starting values for the regressors.
model <- BE_BA_AREA~ELEV+SPR_02_TC1+SPR_02_TC2+SPT_02_TC3

##Note, precision of 0 is same as flat prior alternatively if you want a flat prior
##do not include the beta.prior tag in the beta.control list.
beta.prior <- prior(dist="NORMAL", mu=rep(0, 5), precision=diag(0, 5))

beta.control <- list(beta.update.method="gibbs", beta.starting=rep(0, 5), beta.prior=beta.prior)


ggt.out <- ggt(formula=model, data=BEF.subset,
               var.update.control=var.update.control,
               beta.control=beta.control,
               var.function = var.function,
               run.control=run.control, DIC=TRUE, verbose=TRUE)

print(ggt.out$accept)
print(ggt.out$DIC)

##Get the effective range of spatial dependence.
eff.range <- -log(0.05)/ggt.out$p.samples[,"phi"]

##Summarize and  plot the chain with coda.
mcmc.obj <- mcmc(cbind(ggt.out$p.samples, eff.range))
summary(mcmc.obj)
##plot(mcmc.obj)

##############################################
##ggt setup and model fit with single
##block update of the variance and regression
##parameters
##############################################

var.update.control <- list("sigmasq"=list(sample.order=0, starting=30, tuning=0.5, prior=sigmasq.prior),
                           "tausq"=list(sample.order=0, starting=30, tuning=0.5, prior=tausq.prior),
                           "phi"=list(sample.order=0, starting=0.006, tuning=0.8, prior=phi.prior))

#Using Metropolis-Hastings and default flat prior for regression parameters.
beta.control <- list(beta.update.method="mh", sample.order=0)

ggt.out <- ggt(formula=model, data=BEF.subset,
               var.update.control=var.update.control,
               beta.control=beta.control,
               var.function = var.function,
               run.control=run.control, DIC=TRUE, verbose=TRUE)

##Note single acceptance rate.
print(ggt.out$accept)

}

\keyword{misc}
