\name{spDynLM}
\alias{spDynLM}
\title{Function for fitting univariate Bayesian dynamic space-time regression models}

\description{
  The function \code{spDynLM} fits Gaussian univariate Bayesian
  dynamic space-time regression models for settings where space is viewed as continuous but time is taken
to be discrete.
}

\usage{
spDynLM(formula, data = parent.frame(), coords, knots,
      starting, tuning, priors, cov.model, get.fitted=FALSE, 
      n.samples, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a list of \eqn{N_t}{N_t} symbolic regression models to be fit. Each model represents a time step. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spDynLM} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta}, \code{sigma.sq},
    \code{tau.sq}, \code{phi}, \code{nu}, and \code{sigma.eta}.  The value portion of each tag is the parameter's starting value. }
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing)
    or a vector of length two or three with the first and second elements recording the
    number of columns and rows in the desired knot grid. The third,
    optional, element sets the offset of the outermost knots from the
    extent of the \code{coords}. } 
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{phi} and \code{nu}. The value portion of each tag defines the variance of the Metropolis sampler Normal proposal distribution.}

  \item{priors}{a list with tags \code{beta.0.norm}, \code{sigma.sq.ig},
    \code{tau.sq.ig}, \code{phi.unif}, \code{nu.unif}, and \code{sigma.eta.iw}. Variance parameters, \code{simga.sq} and
    \code{tau.sq}, are assumed to follow an
    inverse-Gamma distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform
    distributions. The \code{beta.0.norm} is a multivariate Normal distribution with
  hyperparameters passed as a list of length two with the first
  and second elements corresponding to the mean vector and positive
  definite covariance matrix, respectively. The hyperparameters of the inverse-Wishart, \code{sigma.eta.iw},  are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{p\times p}{pxp} \emph{scale} matrix,
    respectively. The inverse-Gamma hyperparameters are passed in a list with two vectors that hold the \emph{shape} and \emph{scale}, respectively. The Uniform hyperparameters are passed in a list with two vectors that hold the lower and upper support values, respectively.}
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{get.fitted}{if \code{TRUE}, posterior predicted and fitted
    values are collected. Note, posterior predicted samples are only
    collected for those \eqn{y_t(s)}{y_t(s)} that are \code{NA}.}
    \item{n.samples}{the number of MCMC iterations.}
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report Metropolis sampler acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{spDynLM}, which is a list with the following
  tags:
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
 
  \item{p.theta.samples}{a \code{coda} object of posterior samples for
    \eqn{\tau^2_t}{tau^2_t}, \eqn{\sigma^2_t}{sigma^2_t}, \eqn{\phi_t}{phi_t}, \eqn{\nu_t}{nu_t}. }
  
  \item{p.beta.0.samples}{a \code{coda} object of posterior samples for
    \eqn{\beta}{\beta} at \eqn{t=0}{t=0}.}
  
  \item{p.beta.samples}{a \code{coda} object of posterior samples for
    \eqn{\beta_t}{\beta_t}.}
  
  \item{p.sigma.eta.samples}{a \code{coda} object of posterior samples for
    \eqn{\Sigma_\eta}{\Sigma_\eta}.}
  
  \item{p.u.samples}{a \code{coda} object of posterior samples for
    spatio-temporal random effects \eqn{u}{u}. Samples are over the
    columns and time steps increase in blocks of \eqn{n}{n} down the
    columns, e.g., the first \eqn{n}{n} rows correspond to locations
    \eqn{1,2, \ldots, n}{1, 2, ..., n} in \eqn{t=1}{t=1} and the last \eqn{n}{n} rows correspond to locations
    \eqn{1,2, \ldots, n}{1, 2, ..., n} in \eqn{t=N_t}{t=N_t}.}
  
  \item{p.y.samples}{a \code{coda} object of posterior samples for
    \eqn{y}{y}. If \eqn{y_t(s)}{y_t(s)} is specified as \code{NA} the
    \code{p.y.samples} hold the associated posterior predictive
    samples. Samples are over the
    columns and time steps increase in blocks of \eqn{n}{n} down the
    columns, e.g., the first \eqn{n}{n} rows correspond to locations
    \eqn{1,2, \ldots, n}{1, 2, ..., n} in \eqn{t=1}{t=1} and the last \eqn{n}{n} rows correspond to locations
    \eqn{1,2, \ldots, n}{1, 2, ..., n} in \eqn{t=N_t}{t=N_t}.} 	
  
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\details{
  Suppose, \eqn{y_t(s)}{y_t(s)} denotes the observation at location \eqn{s}{s} and time
\eqn{t}{t}. We model \eqn{y_t(s)}{y_t(s)} through a \emph{measurement equation} that
provides a regression specification with a space-time varying
intercept and  serially and spatially
uncorrelated zero-centered Gaussian disturbances as measurement error
\eqn{\epsilon_t(s)}{epsilon_t(s)}. Next a \emph{transition equation}
  introduces a \eqn{p\times 1}{px1} coefficient vector, say \eqn{\beta_t}{beta_t}, which is a purely
temporal component (i.e., time-varying regression parameters), and a
spatio-temporal component \eqn{u_t(s)}{u_t(s)}. Both these are generated through
transition equations, capturing their Markovian dependence in
time. While the transition equation of the purely temporal component
is akin to usual state-space modeling, the spatio-temporal component is
generated using Gaussian spatial processes. The overall model is written
as
\deqn{y_t(s)  = x_t(s)'\beta_t + u_t(s) + \epsilon_t(s), t=1,2,\ldots,N_t}{y_t(s) =
  x_t(s)'\beta_t + u_t(s) + \epsilon_t(s), t=1,2,..,N_t}

\deqn{\epsilon_t(s) \sim N(0,\tau_{t}^2)}{\epsilon_t(s) ~ N(0,\tau_{t}^2)}

\deqn{\beta_t  = \beta_{t-1} + \eta_t; \eta_t \sim
  N(0,\Sigma_{\eta})}{\beta_t  = \beta_{t-1} + \eta_t; \eta_t ~ N(0,\Sigma_{\eta})}

\deqn{u_t(s) = u_{t-1}(s) + w_t(s); w_t(s) \sim GP(0, C_t(\cdot,\theta_t))}{u_t(s) = u_{t-1}(s) + w_t(s); w_t(s) ~ GP(0, C_t(.,\theta_t))}

Here \eqn{x_t(s)}{x_t(s)} is a \eqn{p\times 1}{px1} vector of predictors
and \eqn{\beta_t}{beta_t} is a \eqn{p\times 1}{px1} vector of
coefficients. In addition to an intercept, \eqn{x_t(s)}{x_t(s)} can
include location specific variables useful for explaining the
variability in \eqn{y_t(s)}{y_t(s)}. The \eqn{GP(0,
  C_t(\cdot,\theta_t))}{GP(0, C_t(.,\theta_t))} denotes a spatial
Gaussian process with covariance function
\eqn{C_{t}(\cdot;\theta_t)}{C_t(.;\theta_t)}. We specify
\eqn{C_{t}(s_1,s_2;\theta_t)=\sigma_t^2\rho(s_1,s_2;\phi_t)}{C_t(s_1,s_2;\theta_t)=\sigma_t^2\rho(s_1,s_2;\phi_t,
\nu_t)}, where \eqn{\theta_t = \{\sigma_t^2,\phi_t,\nu_t\}}{\theta_t =
{\sigma_t^2,\phi_t,\nu_t}} and \eqn{\rho(\cdot;\phi)}{\rho(.;\phi)} is a
\emph{correlation function} with \eqn{\phi}{\phi} controlling the
correlation decay and \eqn{\sigma_t^2}{\sigma_t^2} represents the
spatial variance component. The spatial smoothness parameter,
\eqn{\nu}{\nu}, is used if the Matern spatial correlation function is chosen. We further assume \eqn{\beta_0 \sim N(m_0, \Sigma_0)}{\beta_0 ~ N(m_0, \Sigma_0)} and \eqn{u_0(s) \equiv 0}{u_0(s) = 0}, which completes the prior specifications leading to a well-identified Bayesian hierarhical model and also yield reasonable dependence structures. 
}
  
\seealso{
  \code{\link{spLM}}
}

\references{
  Finley, A.O., S. Banerjee, and A.E. Gelfand. (2012) Bayesian dynamic modeling for large space-time datasets using Gaussian predictive processes. \emph{Journal of Geographical Systems}, 14:29--47.

    Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.
  
  Gelfand, A.E., S. Banerjee, and D. Gamerman (2005) Spatial Process
  Modelling for Univariate and Multivariate Dynamic Spatial Data,
  \emph{Environmetrics}, 16:465--479. 
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
data("NETemp.dat")
ne.temp <- NETemp.dat

set.seed(1)

##take a chunk of New England
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3000000,]

##subset first 2 years (Jan 2000 - Dec. 2002)
y.t <- ne.temp[,4:27]
N.t <- ncol(y.t) ##number of months
n <- nrow(y.t) ##number of observation per months

##add some missing observations to illistrate prediction
miss <- sample(1:N.t, 10)
holdout.station.id <- 5
y.t.holdout <- y.t[holdout.station.id, miss]
y.t[holdout.station.id, miss] <- NA

##scale to km
coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))

##set starting and priors
p <- 2 #number of regression parameters in each month

starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=diag(rep(0.01, p)))

tuning <- list("phi"=rep(5, N.t)) 

priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
               "sigma.eta.IW"=list(2, diag(0.001,p)))

##make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)

n.samples <- 2000

m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
               starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
               cov.model="exponential", n.samples=n.samples, n.report=25) 

burn.in <- floor(0.75*n.samples)

quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}

beta <- apply(m.1$p.beta.samples[burn.in:n.samples,], 2, quant)
beta.0 <- beta[,grep("Intercept", colnames(beta))]
beta.1 <- beta[,grep("elev", colnames(beta))]

plot(m.1$p.beta.0.samples)

par(mfrow=c(2,1))
plot(1:N.t, beta.0[1,], pch=19, cex=0.5, xlab="months", ylab="beta.0", ylim=range(beta.0))
arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[3,], length=0.02, angle=90)
arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[2,], length=0.02, angle=90)

plot(1:N.t, beta.1[1,], pch=19, cex=0.5, xlab="months", ylab="beta.1", ylim=range(beta.1))
arrows(1:N.t, beta.1[1,], 1:N.t, beta.1[3,], length=0.02, angle=90)
arrows(1:N.t, beta.1[1,], 1:N.t, beta.1[2,], length=0.02, angle=90)

theta <- apply(m.1$p.theta.samples[burn.in:n.samples,], 2, quant)
sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
tau.sq <- theta[,grep("tau.sq", colnames(theta))]
phi <- theta[,grep("phi", colnames(theta))]

par(mfrow=c(3,1))
plot(1:N.t, sigma.sq[1,], pch=19, cex=0.5, xlab="months", ylab="sigma.sq", ylim=range(sigma.sq))
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[3,], length=0.02, angle=90)
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[2,], length=0.02, angle=90)

plot(1:N.t, tau.sq[1,], pch=19, cex=0.5, xlab="months", ylab="tau.sq", ylim=range(tau.sq))
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[3,], length=0.02, angle=90)
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[2,], length=0.02, angle=90)

plot(1:N.t, 3/phi[1,], pch=19, cex=0.5, xlab="months", ylab="eff. range (km)", ylim=range(3/phi))
arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[3,], length=0.02, angle=90)
arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[2,], length=0.02, angle=90)

y.hat <- apply(m.1$p.y.samples[,burn.in:n.samples], 1, quant)
y.hat.med <- matrix(y.hat[1,], ncol=N.t)
y.hat.up <- matrix(y.hat[3,], ncol=N.t)
y.hat.low <- matrix(y.hat[2,], ncol=N.t)

y.obs <- as.vector(as.matrix(y.t[-holdout.station.id, -miss]))
y.obs.hat.med <- as.vector(y.hat.med[-holdout.station.id, -miss])
y.obs.hat.up <- as.vector(y.hat.up[-holdout.station.id, -miss])
y.obs.hat.low <- as.vector(y.hat.low[-holdout.station.id, -miss])

y.ho <- as.matrix(y.t.holdout)
y.ho.hat.med <- as.vector(y.hat.med[holdout.station.id, miss])
y.ho.hat.up <- as.vector(y.hat.up[holdout.station.id, miss])
y.ho.hat.low <- as.vector(y.hat.low[holdout.station.id, miss])

par(mfrow=c(2,1))
plot(y.obs, y.obs.hat.med, pch=19, cex=0.5, xlab="observed",
ylab="fitted", main="Observed vs. fitted")
arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.up, length=0.02, angle=90)
arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.low, length=0.02, angle=90)
lines(-50:50, -50:50, col="blue")

plot(y.ho, y.ho.hat.med, pch=19, cex=0.5, xlab="observed",
ylab="predicted", main="Observed vs. predicted")
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.up, length=0.02, angle=90)
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.low, length=0.02, angle=90)
lines(-50:50, -50:50, col="blue")
}
}

\keyword{model}