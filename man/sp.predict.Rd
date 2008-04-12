\name{sp.predict}
\alias{sp.predict}          
\title{Prediction for new points given a ggt.sp, sp.lm, or bayes.geostat.exact object}

\description{
  The function \code{sp.predict} collects a posterior predictive sample
  for a set of new points given a \code{\link{ggt.sp}},
  \code{\link{sp.lm}}, or \code{\link{bayes.geostat.exact}} object.
}

\usage{
  sp.predict(sp.obj, pred.coords, pred.covars,
             start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{ggt.sp}}, \code{\link{bayes.geostat.exact}}, or \code{\link{sp.lm}}}
  \item{pred.coords}{an \eqn{n \times 2}{n x 2} matrix of \eqn{m} prediction
    point coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
    The first column is assumed to be easting coordinates and the second column
    northing coordinates. }
  \item{pred.covars}{an \eqn{n \times p}{n x p} design matrix associated
    with the new points.  If this is a multivariate prediction defined
    by \eqn{m}{m} models, the multivariate design matrix can be created
    by passing a list of the \eqn{m}{m} univariate design matrices to
    the \code{\link{mk.mv.X}} function. }
  \item{start}{specifies the first sample included in the prediction calculation.
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
  \item{obs.coords}{the matrix of the observation coordinates.  }
  \item{pred.coords}{the matrix of prediction point coordinates specified by \code{pred.coords}. }
  \item{pp.samples}{a matrix that holds samples from the posterior
    predictive distribution(s). For \code{ggt.sp}, the rows of this matrix
    correspond to the predicted points and the columns are the posterior
    predictive samples.  If prediction is for \eqn{m}{m} response variables the \code{pp.samples} matrix has
    \eqn{mn}{mn} rows.  The predictions for points are held in rows
    1:m, (m+1):2m, ..., ((i-1)m+1):im, ..., ((n-1)m+1):nm, where i
    = 1 ... n (e.g., the samples for the first point are in rows 1:m,
    second point in rows (m+1):2m, etc.). For
    \code{bayes.geostat.exact}, the rows of this matrix
    correspond to the predicted points and the columns are the posterior
    predictive samples.}
}

\details{
  Please refer to Section 3.2 in the vignette.
}

\seealso{\code{\link{ggt.sp}}, \code{\link{bayes.geostat.exact}}, \code{\link{sp.lm}}}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

  Banerjee, S., Gelfand, A.E., Finley, A.O., and Sang, H. (In press). Gaussian predictive process models for large spatial datasets. \emph{Journal of the Royal Statistical Society Series B}.
  %%Further information on the package \pkg{spBayes} can be found at: \url{http://blue.fr.umn.edu/spatialBayes}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}.
}

\examples{
\dontrun{
##Portions of this example requires MBA package to make surfaces
library(MBA)

###########################################
##  Prediction for bayes.geostat.exact
###########################################
data(FBC07.dat)
Y <- FBC07.dat[1:150,"Y.2"]
coords <- as.matrix(FBC07.dat[1:150,c("coord.X", "coord.Y")])

n.samples <-1000
n = length(Y)
p = 1

phi <- 0.15
nu <- 0.5

beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

alpha <- 5/5

sigma.sq.prior.shape <- 2.0
sigma.sq.prior.rate <- 5.0

##############################
##Simple linear model with
##the default exponential
##spatial decay function
##############################
m.1 <- bayes.geostat.exact(Y~1, n.samples=n.samples,
                           beta.prior.mean=beta.prior.mean,
                           beta.prior.precision=beta.prior.precision,
                           coords=coords, phi=phi, alpha=alpha,
                           sigma.sq.prior.shape=sigma.sq.prior.shape,
                           sigma.sq.prior.rate=sigma.sq.prior.rate,
                           sp.effects=TRUE)

##Now prediction
set.seed(1)
pred.coords <- expand.grid(seq(0,100,length=10),seq(0,100,length=10))
pred.covars <- as.matrix(rep(1,nrow(pred.coords)))

m.1.pred <- sp.predict(m.1, pred.coords=pred.coords,
                       pred.covars=pred.covars, thin=5)

par(mfrow=c(2,2))
obs.surf <-
  mba.surf(cbind(coords, Y), no.X=100, no.Y=100, extend=T)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords, pch=19, cex=1, col="green")
contour(obs.surf, add=T)

w.hat <- rowMeans(m.1$sp.effects)
w.surf <-
  mba.surf(cbind(coords, w.hat), no.X=100, no.Y=100, extend=T)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Random effects")
points(coords, pch=19, cex=1, col="green")
contour(w.surf, add=T)

y.hat <- rowMeans(m.1.pred)
y.surf <-
  mba.surf(cbind(pred.coords, y.hat), no.X=100, no.Y=100, extend=T)$xyz.est
image(y.surf, xaxs = "r", yaxs = "r", main="Predicted response")
points(pred.coords, pch=19, cex=1, col="black")
rect(0, 0, 50, 50, col=NA, border="green")
contour(y.surf, add=T)

y.var <- apply(m.1.pred, 1, var)
y.surf <-
  mba.surf(cbind(pred.coords, y.var), no.X=100, no.Y=100, extend=T)$xyz.est
image(y.surf, xaxs = "r", yaxs = "r", main="Predicted response\nvariance")
points(coords, pch=19, cex=1, col="green")
points(pred.coords, pch=19, cex=1, col="black")
rect(0, 0, 50, 50, col=NA, border="green")
contour(y.surf, add=T)
 
###########################################
##       Prediction for sp.lm
###########################################
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
w <- rf.n200.dat$w

pred.coords <- expand.grid(seq(1,10,1), seq(1,10,1))
n.pred <- nrow(pred.coords)

###############################
##Prediction with a sp.lm model
###############################
m.2 <- sp.lm(Y~1, coords=coords,
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100, sp.effects=TRUE)

pred <- sp.predict(m.2, pred.coords,
                   pred.covars=as.matrix(rep(1,n.pred)))

par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, Y), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

y.hat <- rowMeans(pred$y.pred)
y.pred.surf <-
  mba.surf(cbind(pred.coords, y.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(y.pred.surf, xaxs = "r", yaxs = "r", main="Predicted response")
points(coords, pch=1, cex=1)
points(pred.coords, pch=19, cex=1)
contour(y.pred.surf, add=T)
legend(1.5,2.5, legend=c("Obs.", "Pred."), pch=c(1,19),
       cex=c(1,1), bg="white")

###############################
##Prediction with a sp.lm
##predictive process model
###############################
m.3 <- sp.lm(Y~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=2000, verbose=TRUE, n.report=100, sp.effects=TRUE)

print(summary(m.3$p.samples))
plot(m.3$p.samples)

pred <- sp.predict(m.3, pred.coords,
                   pred.covars=as.matrix(rep(1,n.pred)))

par(mfrow=c(1,2))
obs.surf <-
  mba.surf(cbind(coords, Y), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)

y.hat <- rowMeans(pred$y.pred)
y.pred.surf <-
  mba.surf(cbind(pred.coords, y.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(y.pred.surf, xaxs = "r", yaxs = "r", main="Predicted response")
points(coords, pch=1, cex=1)
points(m.3$knot.coords, pch=3, cex=1)
points(pred.coords, pch=19, cex=1)
contour(y.pred.surf, add=T)
legend(1.5,2.5, legend=c("Obs.", "Knots", "Pred."),
       pch=c(1,3,19), cex=c(1,1,1), bg="white")

###########################################
##       Prediction for ggt.sp
###########################################
data(FBC07.dat)

##Divide the data into model and prediction sets
Y.1 <- FBC07.dat[1:100,"Y.1"]
Y.2 <- FBC07.dat[1:100,"Y.2"]
model.coords <- as.matrix(FBC07.dat[1:100,c("coord.X", "coord.Y")])
pred.coords <- as.matrix(FBC07.dat[151:200,c("coord.X", "coord.Y")])

#############################
##   Univariate model
#############################

##Fit some model with ggt.sp.
K.prior <- prior(dist="IG", shape=2, scale=5)
Psi.prior <- prior(dist="IG", shape=2, scale=5)
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.5, prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=1000)

Fit <- ggt.sp(formula=Y.2~1, run.control=run.control,
              coords=model.coords,
              var.update.control=var.update.control,
              beta.update.control=beta.control,
              cov.model="exponential")

##Now make predictions for the holdout set.
##Step 1. make the design matrix for the prediction points.
pred.covars <- as.matrix(rep(1, nrow(pred.coords)))

##Step 2. call sp.predict.
Pred <- sp.predict(Fit, pred.covars=pred.covars,
                   pred.coords=pred.coords)

##Step 3. check out the predicted random effects and
##predicted response variable.

Pred.sp.effects.surf <-
  mba.surf(cbind(pred.coords, rowMeans(Pred$pred.sp.effects)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

Pred.Y.surf <-
  mba.surf(cbind(pred.coords, rowMeans(Pred$pred.y)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

par(mfrow=c(1,2))
image(Pred.sp.effects.surf, xaxs="r", yaxs="r",
      main="Predicted random spatial effects")
contour(Pred.sp.effects.surf, add=TRUE)

image(Pred.Y.surf, xaxs="r", yaxs="r",
      main="Predicted Y.2")
contour(Pred.Y.surf, add=TRUE)

#############################
##   Multivariate models
#############################

##Fit some model with ggt.sp.
K.prior <- prior(dist="IWISH", df=2, S=diag(c(3, 6)))
Psi.prior <- prior(dist="IWISH", df=2, S=diag(c(7, 5)))
phi.prior <- prior(dist="UNIF", a=0.06, b=3)

K.starting <- matrix(c(2,-3, 0, 1), 2, 2)
Psi.starting <- diag(c(3, 2))

var.update.control <-
  list("K"=list(starting=K.starting, tuning=diag(c(0.1, 0.5, 0.1)),
         prior=K.prior),
       "Psi"=list(starting=Psi.starting, tuning=diag(c(0.1, 0.5, 0.1)),
         prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5,
         prior=list(phi.prior, phi.prior))
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=1000, "sp.effects"=FALSE)

Fit.mv <-
  ggt.sp(formula=list(Y.1~1, Y.2~1), run.control=run.control,
         coords=model.coords,
         var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

##Now make predictions for the holdout set.
##Step 1. make the design matrix for the prediction points using
##the mk.mv.X function.
pred.covars.1 <- as.matrix(rep(1, nrow(pred.coords)))
pred.covars.2 <- as.matrix(rep(1, nrow(pred.coords)))

pred.covars.mv <- mk.mv.X(list(pred.covars.1, pred.covars.2))

##Step 2. call sp.predict.
Pred.mv <- sp.predict(Fit.mv, pred.covars=pred.covars.mv,
                      pred.coords=pred.coords)

##Step 3. check out the predicted random effects and
##predicted response variables.   Recall, these are
##organized as m consecutive rows for each point.
Pred.sp.effects.1 <-
  Pred.mv$pred.sp.effects[seq(1, nrow(Pred.mv$pred.sp.effects), 2),]

Pred.sp.effects.2 <-
  Pred.mv$pred.sp.effects[seq(2, nrow(Pred.mv$pred.sp.effects), 2),]

Pred.Y.1 <-
  Pred.mv$pred.sp.effects[seq(1, nrow(Pred.mv$pred.y), 2),]

Pred.Y.2 <-
  Pred.mv$pred.sp.effects[seq(2, nrow(Pred.mv$pred.y), 2),]

##Then map.
Pred.sp.effects.1.surf <-
  mba.surf(cbind(pred.coords, rowMeans(Pred.sp.effects.1)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

Pred.sp.effects.2.surf <-
  mba.surf(cbind(pred.coords, rowMeans(Pred.sp.effects.2)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

Pred.Y.1.surf <-
  mba.surf(cbind(pred.coords, rowMeans(Pred.Y.1)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

Pred.Y.2.surf <-
  mba.surf(cbind(pred.coords, rowMeans(Pred.Y.2)),
           no.X=100, no.Y=100, extend=TRUE)$xyz.est

par(mfrow=c(2,2))
image(Pred.sp.effects.surf, xaxs="r", yaxs="r",
      main="Predicted random spatial effects Y.1")
contour(Pred.sp.effects.1.surf, add=TRUE)

image(Pred.sp.effects.surf, xaxs="r", yaxs="r",
      main="Predicted random spatial effects Y.2")
contour(Pred.sp.effects.2.surf, add=TRUE)

image(Pred.sp.effects.surf, xaxs="r", yaxs="r",
      main="Predicted Y.1")
contour(Pred.Y.1.surf, add=TRUE)

image(Pred.sp.effects.surf, xaxs="r", yaxs="r",
      main="Predicted Y.2")
contour(Pred.Y.2.surf, add=TRUE)

}
}

\keyword{misc}
