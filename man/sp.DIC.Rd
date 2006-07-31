\name{sp.DIC}
\alias{sp.DIC}          
\title{Calculates model DIC and associated statistics given a ggt.sp object}

\description{
  The function \code{sp.DIC} calculates model DIC and associated statistics given a \code{\link{ggt.sp}} object.
}

\usage{
  sp.DIC(ggt.sp.obj, DIC.marg=TRUE, DIC.unmarg=TRUE,
         start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{ggt.sp.obj}{an object returned by \code{\link{ggt.sp}} (i.e.,
    object of class \code{ggt.sp}). }
  \item{DIC.marg}{a logical value indicating if marginalized DIC and
    associated statistics should be calculated. }
  \item{DIC.unmarg}{a logical value indicating if unmarginalized DIC and
    associated statistics should be calculated. }
  \item{start}{specifies the first sample included in the DIC calculation.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{end}{specifies the last sample included  in the DIC calculation.
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
  \item{DIC.marg}{a matrix holding marginalized DIC and associated
    statistics.  }
  \item{DIC.unmarg}{a matrix holding unmarginalized DIC and associated
    statistics.  }
  \item{sp.effects}{if DIC.ummarg is true and if the \code{ggt.sp.obj}
    does not include \code{sp.effects} then \code{sp.DIC} calculates the
    random spatial effects and includes them in the return object. }
}

\details{
  Please refer to Section 3.3 in the vignette.
}

\seealso{\code{\link{ggt.sp}}}

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
data(FBC07)

Y.2 <- FBC07[1:100,"Y.2"]
coords <- as.matrix(FBC07[1:100,c("coord.X", "coord.Y")])

##Fit some model with ggt.sp.
K.prior <- prior(dist="IG", shape=2, scale=5)
Psi.prior <- prior(dist="IG", shape=2, scale=5)
phi.prior <- prior(dist="LOGUNIF", a=0.06, b=3)

var.update.control <-
  list("K"=list(starting=5, tuning=0.5, prior=K.prior),
       "Psi"=list(starting=5, tuning=0.5, prior=Psi.prior),
       "phi"=list(starting=0.1, tuning=0.5, prior=phi.prior)
       )

beta.control <- list(update="GIBBS", prior=prior(dist="FLAT"))

run.control <- list("n.samples"=1000, "sp.effects"=TRUE)

Fit <-
  ggt.sp(formula=Y.2~1, run.control=run.control,
         coords=coords, var.update.control=var.update.control,
         beta.update.control=beta.control,
         cov.model="exponential")

##Now with the ggt.sp object, Fit, calculate the DIC
##for both the unmarginalized and marginalized models.
##The likelihoods for these models are given by equation 6 and 7
##within the vignette.

DIC <- sp.DIC(Fit)
print(DIC)

}

\keyword{misc}
