\name{spPredict}
\alias{spPredict}          
\title{Prediction for new points given a spGGT, spLM, spMvLM, spGLM,
  spMvGLM, or bayesGeostatExact object}

\description{
  The function \code{spPredict} collects a posterior predictive sample
  for a set of new points given a \code{\link{spGGT}},
  \code{\link{spLM}}, \code{\link{spMvLM}}, \code{\link{spGLM}},
  \code{\link{spMvGLM}}, or \code{\link{bayesGeostatExact}} object.
}

\usage{
  spPredict(sp.obj, pred.coords, pred.covars,
             start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{spGGT}},
  \code{\link{bayesGeostatExact}}, \code{\link{spLM}},
  \code{\link{spMvLM}}, \code{\link{spGLM}}, or \code{\link{spMvGLM}}}
  \item{pred.coords}{an \eqn{n \times 2}{n x 2} matrix of \eqn{m} prediction
    point coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
    The first column is assumed to be easting coordinates and the second column
    northing coordinates. }
  \item{pred.covars}{an \eqn{n \times p}{n x p} design matrix associated
    with the new points.  If this is a multivariate prediction defined
    by \eqn{m}{m} models, the multivariate design matrix can be created
    by passing a list of the \eqn{m}{m} univariate design matrices to
    the \code{\link{mkMvX}} function. }
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
    predictive distribution(s). For multivariate models the rows of this matrix
    correspond to the predicted points and the columns are the posterior
    predictive samples.  If prediction is for \eqn{m}{m} response variables the \code{pp.samples} matrix has
    \eqn{mn}{mn} rows.  The predictions for points are held in rows
    1:m, (m+1):2m, ..., ((i-1)m+1):im, ..., ((n-1)m+1):nm, where i
    = 1 ... n (e.g., the samples for the first point are in rows 1:m,
    second point in rows (m+1):2m, etc.). For
    \code{bayesGeostatExact}, the rows of this matrix
    correspond to the predicted points and the columns are the posterior
    predictive samples.}
}

\seealso{\code{\link{spGGT}}, \code{\link{bayesGeostatExact}},
  \code{\link{spLM}}, \code{\link{spMvLM}}, \code{\link{spGLM}}, \code{\link{spMvGLM}}}

\references{
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.
  
  Finley, A.O,. H. Sang, S. Banerjee, and A.E. Gelfand. (2008). Improving the performance of predictive process modeling for
  large datasets. \emph{Computational Statistics and Data Analysis}, DOI: 10.1016/j.csda.2008.09.008    
    
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

  Banerjee, S., Gelfand, A.E., Finley, A.O., and Sang, H. (In press). Gaussian predictive process models for large spatial datasets. \emph{Journal of the Royal Statistical Society Series B}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}.
}

\examples{
\dontrun{
##Portions of this example requires MBA package to make surfaces
library(MBA)

#########################
##Prediction for spMvGLM
#########################

##Generate some count data
n <- 100
q <- 3   
nltr <- q*(q-1)/2+q

coords <- cbind(runif(n,1,100),runif(n,1,100))

theta <- rep(3/50,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rnorm(nltr,1,1)
K <- A\%*\%t(A)

Psi <- diag(0,q)

c1 <- mvCovInvLogDet(coords=coords, knots=coords,
                     cov.model="exponential",
                     V=K, Psi=Psi, theta=theta,
                     modified.pp=TRUE, SWM=FALSE)

w <- mvrnorm(1,rep(0,nrow(c1$C)),c1$C)

X <- mkMvX(list(matrix(1,n,1), matrix(1,n,1), matrix(1,n,1)))
beta <- c(-1,0,1)
y <- rpois(n*q, exp(X\%*\%beta+w))

y.1 <- y[seq(1,length(y),q)]
y.2 <- y[seq(2,length(y),q)]
y.3 <- y[seq(3,length(y),q)]

##Specify starting values and collect samples. For
##a true analysis, several longer chains should be
##run.
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

beta.starting <- coefficients(glm(y~X-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~X-1, family="poisson"))))

n.samples <- 5000

m.1 <- spMvGLM(list(y.1~1,y.2~1,y.3~1), family="poisson", coords=coords,
               knots=c(8,8,0),
               starting=
               list("beta"=beta.starting, "phi"=rep(0.06,q),
                    "A"=A.starting, "w"=0), 
               tuning=
               list("beta"=beta.tuning, "phi"=rep(0.01,q),
                    "A"=rep(0.005,nltr), "w"=0.001),
               priors=
               list("beta.Flat", "phi.Unif"=rep(c(0.03, 0.3),q),
                    "K.IW"=list(q+1, diag(0.1,q))),
               cov.model="exponential",
               n.samples=n.samples, verbose=TRUE, n.report=500)


pred.coords <- expand.grid(seq(0,100,10), seq(0,100,10))
m <- nrow(pred.coords)
pred.covars <- mkMvX(list(matrix(1,m,1), matrix(1,m,1), matrix(1,m,1)))

out <-
  spPredict(m.1, start=2000, pred.coords=pred.coords, pred.covars=pred.covars)

y.hat <- rowMeans(out$y.pred)

y.hat.1 <- y.hat[seq(1,length(y.hat),q)]
y.hat.2 <- y.hat[seq(2,length(y.hat),q)]
y.hat.3 <- y.hat[seq(3,length(y.hat),q)]

##Take a look
par(mfrow=c(3,2))
surf <- mba.surf(cbind(coords,y.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.1, cex=1, col="blue")

surf <- mba.surf(cbind(pred.coords,y.hat.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, add=TRUE)
contour(surf, drawlabels=FALSE, add=TRUE)
text(pred.coords, labels=round(y.hat.1,0), cex=1, col="green")

surf <- mba.surf(cbind(coords,y.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.2, cex=1, col="blue")

surf <- mba.surf(cbind(pred.coords,y.hat.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(pred.coords, labels=round(y.hat.2,0), cex=1, col="green")

surf <- mba.surf(cbind(coords,y.3),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.3, cex=1, col="blue")

surf <- mba.surf(cbind(pred.coords,y.hat.3),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(pred.coords, labels=round(y.hat.3,0), cex=1, col="green")


#########################
##Prediction for spGLM
#########################

##Generate count data
n <- 300

coords <- cbind(runif(n,1,100),runif(n,1,100))

phi <- 3/75
sigma.sq <- 3

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- mvrnorm(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
y <- rpois(n, exp(x\%*\%beta+w))

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="poisson"))))
            
n.samples <- 15000

m.1 <- spGLM(y~1, family="poisson", coords=coords, knots=c(8,8,0),
             starting=
             list("beta"=beta.starting, "phi"=0.06,
                  "sigma.sq"=1, "w"=0),
             tuning=
             list("beta"=beta.tuning, "phi"=0.1,
                  "sigma.sq"=0.1, "w"=0.001),
             priors=
             list("beta.Flat", "phi.Unif"=c(0.03, 0.3),
                  "sigma.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=n.samples, verbose=TRUE, n.report=500)

##Make prediction grid
pred.coords <- expand.grid(seq(0,100,5), seq(0,100,5))
out <- spPredict(m.1, start=10000, thin=2, pred.coords=pred.coords,
                 pred.covars=as.matrix(rep(1,nrow(pred.coords))))

y.pred <- rowMeans(out$y.pred)

##Take a look
par(mfrow=c(1,2))
surf <-
  mba.surf(cbind(coords,y),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, add=TRUE)

surf <-
  mba.surf(cbind(pred.coords, y.pred),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Predicted counts")
contour(surf, add=TRUE)
points(pred.coords, pch=19, cex=0.5, col="blue")


#########################
##Prediction for spMvLM
#########################

##Generate some data
n <- 100 ##observed
m <- 50  ##predict
N <- n+m
q <- 3   
nltr <- q*(q-1)/2+q

coords <- cbind(runif(N),runif(N))
theta <- rep(3/0.5,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rnorm(nltr, 5, 1)
V <- A\%*\%t(A)

Psi <- diag(1,q)

c1 <- mvCovInvLogDet(coords=coords, knots=coords,
                     cov.model="exponential",
                     V=V, Psi=Psi, theta=theta)

w <- mvrnorm(1,rep(0,nrow(c1$C)),c1$C)
obs.w <- w[1:(n*q)]

w.1 <- obs.w[seq(1,length(obs.w),q)]
w.2 <- obs.w[seq(2,length(obs.w),q)]
w.3 <- obs.w[seq(3,length(obs.w),q)]

##Specify starting values and collect samples
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
L.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

n.samples <- 1000

obs.coords <- coords[1:n,]

m.1 <- spMvLM(list(w.1~1,w.2~1,w.3~1), coords=obs.coords,
              knots=c(8,8,0),
              starting=list("beta"=rep(1,q), "phi"=rep(3/0.5,q),
                "nu"=rep(1,q), "A"=A.starting, "L"=L.starting), 
              sp.tuning=list("phi"=rep(0.1,q), "nu"=rep(0.1,q),
                "A"=rep(0.1,nltr), "L"=rep(0.1,nltr)),
              priors=list("phi.Unif"=rep(c(3/1,3/0.1),q),
                "nu.Unif"=rep(c(0.1,2),q),
                "K.IW"=list(q+1, diag(1,q)), "Psi.IW"=list(q+1, diag(1,q))),
              modified.pp=FALSE, cov.model="matern",
              n.samples=n.samples, verbose=TRUE, n.report=100)

##Predict for holdout set
pred.coords <- coords[(n+1):nrow(coords),]
pred.covars <- mkMvX(list(matrix(1,m,1), matrix(1,m,1), matrix(1,m,1)))

pred <- spPredict(m.1, pred.coords, pred.covars)

ho.w <- w[(n*q+1):length(w)]
ho.w.1 <- ho.w[seq(1,length(ho.w),q)]
ho.w.2 <- ho.w[seq(2,length(ho.w),q)]
ho.w.3 <- ho.w[seq(3,length(ho.w),q)]

burn.in <- 500

pred.w <- rowMeans(pred$y.pred[,burn.in:ncol(pred$y.pred)])
pred.w.1 <- pred.w[seq(1,length(pred.w),q)]
pred.w.2 <- pred.w[seq(2,length(pred.w),q)]
pred.w.3 <- pred.w[seq(3,length(pred.w),q)]

##Take a look
par(mfrow=c(3,2))
surf <- mba.surf(cbind(obs.coords,w.1),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Observed"); contour(surf, add=TRUE)

surf <- mba.surf(cbind(pred.coords,pred.w.1),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf, main="Predicted"); contour(surf, add=TRUE)
points(m.1$knot.coords, pch=19, cex=1)

surf <- mba.surf(cbind(obs.coords,w.2),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf); contour(surf, add=TRUE)

surf <- mba.surf(cbind(pred.coords,pred.w.2),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf); contour(surf, add=TRUE)
points(m.1$knot.coords, pch=19, cex=1)

surf <- mba.surf(cbind(obs.coords,w.3),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf); contour(surf, add=TRUE)

surf <- mba.surf(cbind(pred.coords,pred.w.3),
                 no.X=100, no.Y=100, extend=T)$xyz.est
image(surf); contour(surf, add=TRUE)
points(m.1$knot.coords, pch=19, cex=1)

###########################################
##  Prediction for bayesGeostatExact
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
m.1 <- bayesGeostatExact(Y~1, n.samples=n.samples,
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

m.1.pred <- spPredict(m.1, pred.coords=pred.coords,
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
##       Prediction for spLM
###########################################
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
w <- rf.n200.dat$w

pred.coords <- expand.grid(seq(1,10,1), seq(1,10,1))
n.pred <- nrow(pred.coords)

###############################
##Prediction with a spLM model
###############################
m.2 <- spLM(Y~1, coords=coords,
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100)

pred <- spPredict(m.2, pred.coords,
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
legend(1.5,2.5, legend=c("Obs.", "Pred."), pch=c(1,19), bg="white")

###############################
##Prediction with a spLM
##predictive process model
###############################
m.3 <- spLM(Y~1, coords=coords, knots=c(6,6,0),
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=2000, verbose=TRUE, n.report=100)

print(summary(m.3$p.samples))
plot(m.3$p.samples)

pred <- spPredict(m.3, pred.coords,
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
       pch=c(1,3,19), bg="white")

###########################################
##       Prediction for spGGT
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

##Fit some model with spGGT.
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

Fit <- spGGT(formula=Y.2~1, run.control=run.control,
              coords=model.coords,
              var.update.control=var.update.control,
              beta.update.control=beta.control,
              cov.model="exponential")

##Now make predictions for the holdout set.
##Step 1. make the design matrix for the prediction points.
pred.covars <- as.matrix(rep(1, nrow(pred.coords)))

##Step 2. call spPredict.
Pred <- spPredict(Fit, pred.covars=pred.covars,
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

##Fit some model with spGGT.
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
  spGGT(formula=list(Y.1~1, Y.2~1), run.control=run.control,
         coords=model.coords,
         var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

##Now make predictions for the holdout set.
##Step 1. make the design matrix for the prediction points using
##the mkMvX function.
pred.covars.1 <- as.matrix(rep(1, nrow(pred.coords)))
pred.covars.2 <- as.matrix(rep(1, nrow(pred.coords)))

pred.covars.mv <- mkMvX(list(pred.covars.1, pred.covars.2))

##Step 2. call spPredict.
Pred.mv <- spPredict(Fit.mv, pred.covars=pred.covars.mv,
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
