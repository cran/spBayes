\name{sp.predict}
\alias{sp.predict}          
\title{Prediction for new points given a ggt.sp object}

\description{
  The function \code{sp.predict} collects a posterior predictive sample for a set of new points given a \code{\link{ggt.sp}} object.
}

\usage{
  sp.predict(ggt.sp.obj, pred.joint = TRUE, pred.coords, pred.covars, use.covar.names = TRUE, start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{ggt.sp.obj}{an object returned by \code{\link{ggt.sp}} (i.e.,
    object of class \code{ggt.sp}). }
  \item{pred.joint}{if \code{TRUE} the spatial dependence among the
    prediction points is preserved; otherwise, predicted points are
    considered independent.  See below for details. }
  \item{pred.coords}{an \eqn{m \times 2}{m x 2} matrix of \eqn{m} prediction
    point coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
    The first column is assumed to be easting coordinates and the second column
    northing coordinates (i.e., x,y coordinates). }
  \item{pred.covars}{an \eqn{m \times q}{m x q} matrix or data frame containing the
    regressors associated with \code{pred.coords}.  Several scenarios
    exist which determine the specification of this matrix or data frame.  If the model
    contained in \code{ggt.sp.obj} has only an intercept, then
    \code{pred.covars} is optional.  Specifically, in this case, if
    \code{pred.covars} is not specified then \code{ggt.sp.obj$X} is
    checked for only an intercept, if true an intercept is added to
    the \eqn{m} prediction points; otherwise, and error is
    thrown.   If \code{use.covar.names} is \code{TRUE}, then both then column names for \code{ggt.sp.obj$X} and
    \code{pred.covars} are matched.  If in this case an intercept is
    used in \code{ggt.sp.obj$X}, but not specified in \code{pred.covars},
    it is added with a warning.  If one ore more regressor column names are missing
    from \code{pred.covars} an error is thrown.  If \code{use.covar.names} is
    \code{FALSE}, then the column order is assumed to match the intercept
    and regressor columns in \code{ggt.sp.obj$X}.  In this case, it is
    the users responsibility to add an intercept to \code{pred.covars}
    if needed.  Further, an error is thrown if the number of columns in
    \code{pred.covars} is not equal to the number of columns in
    \code{ggt.sp.obj$X}. }
  \item{use.covar.names}{if \code{TRUE} an attempt is made to match
    regressor names between \code{ggt.sp.obj$X} and \code{pred.covars};
    otherwise, column order must match between
    \code{ggt.sp.obj$X} and \code{pred.covars}. }
  \item{start}{specifies the first sample included in the prediction calculation.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{end}{specifies the last sample included  in the prediction calculation.
    The default is to include from \code{start} to
    \code{nrow(ggt.sp.obj$p.samples)}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10}
    then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{verbose}{if \code{TRUE} calculation progress is printed to the
    screen; otherwise, nothing is printed to the screen.  }
  \item{...}{currently no additional arguments.  }
}

\value{
  \item{obs.coords}{the \eqn{n \times 2}{n x 2} matrix of the observation
    coordinates from \code{ggt.sp.obj$coords}.  }
  \item{pred.coords}{the \eqn{m \times 2}{m x 2} matrix of prediction
    point coordinates specified by \code{pred.coords}. }
  \item{pp.samples}{an \eqn{n \times m}{n x m} matrix that holds \eqn{m}
    samples from the posterior predictive distribution of the response
    variable specified in the \code{ggt.sp.obj}. }
}

\seealso{\code{\link{ggt.sp}}, \code{\link{sp.effect}}}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
  
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

#specify the number of samples
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

###############################################
##prediction
###############################################

##resample
n.pred <- 100
BEF.pred <- BEF[sample(1:nrow(BEF),n.pred),]

pred.coords <- cbind(BEF.pred$XUTM, BEF.pred$YUTM)

pred <- sp.predict(ggt.sp.out, pred.joint=TRUE, pred.coords=pred.coords, pred.covars=BEF.pred, start=500, thin=2)

par(mfrow=c(1,2))
int.obs <- interp(BEF.subset$XUTM, BEF.subset$YUTM, BEF.subset$BE_BA_AREA)
image(int.obs, xlab="Longitude (meters)", ylab="Latitude (meters)", main="Observed density\nof American beech")
contour(int.obs, add=TRUE)

int.pred <- interp(BEF.pred$XUTM, BEF.pred$YUTM, rowMeans(pred$pp.samples))
image(int.pred, xlab="Longitude (meters)", ylab="Latitude (meters)", main="Predicted density\nof American beech")
contour(int.pred, add=TRUE)

}

\keyword{misc}
