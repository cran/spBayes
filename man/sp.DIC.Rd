\name{sp.DIC}
\alias{sp.DIC}
\title{Calculates Deviance Information Criterion (DIC) for ggt.sp object}

\description{
  The function \code{sp.DIC} calculates DIC and associated statistics for a given \code{\link{ggt.sp}} object.
}

\usage{
sp.DIC(ggt.sp.obj, start=1, end, thin=1, verbose=TRUE, ...)
}

\arguments{
  \item{ggt.sp.obj}{an object returned by \code{\link{ggt.sp}} (i.e.,
    object of class \code{ggt.sp}). }
  \item{start}{specifies the first sample included in the DIC calculation.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{end}{specifies the last sample included in the DIC calculation.
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
  \item{DIC}{a matrix that holds \eqn{DIC} and associated values \eqn{\bar{D}}{bar(D)},
    \eqn{D[\bar{\theta}]}{D[bar(theta)]}, and \eqn{pD} (Banerjee et al. 2004
    and Spiegelhalter et al. 2002). }
}

\seealso{\code{\link{ggt.sp}}}

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

\keyword{misc}
