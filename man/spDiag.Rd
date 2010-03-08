\name{spDiag}
\alias{spDiag}          
\title{Model fit diagnostics DIC and GP}

\description{
  The function \code{spDiag} calculates DIC, GP, and associated
  statistics given a \code{\link{bayesLMRef}}, \code{\link{bayesLMConjugate}}, \code{\link{spGGT}}, \code{\link{spLM}},  
  \code{\link{spMvLM}}, \code{\link{bayesGeostatExact}},
  \code{\link{spGLM}}, or \code{\link{mvLM}} object.
}

\usage{
  spDiag(sp.obj, start=1, end, thin=1, n.report=100, verbose=TRUE, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{bayesLMRef}}, \code{\link{bayesLMConjugate}}, \code{\link{spGGT}},
    \code{\link{spLM}}, \code{\link{spMvLM}}, \code{\link{mvLM}},
  \code{\link{bayesGeostatExact}}, or \code{\link{spGLM}}}
  \item{start}{specifies the first sample included in the calculation.
    This is useful for those who choose to acknowledge chain burn-in.  }
  \item{end}{specifies the last sample included  in the prediction calculation.
    The default is to use all posterior samples in \code{sp.obj}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10}
    then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{verbose}{if \code{TRUE} calculation progress is printed to the
    screen; otherwise, nothing is printed to the screen.  }
  \item{n.report}{the interval to report progress. } 
  \item{...}{currently no additional arguments.  }
}

\value{
  A list with some of the following tags:
  \item{DIC}{a matrix holding DIC and associated statistics, see
  Banerjee et al. (2004) for details.}
  \item{GP}{a matrix holding GP and associated statistics, see Gelfand
  and Ghosh (1998) for details. This is only
  available for \code{\link{bayesLMRef}},
  \code{\link{bayesLMConjugate}},\code{\link{spLM}},
  \code{\link{spMvLM}},
  \code{\link{bayesGeostatExact}}, and \code{\link{mvLM}} objects}
  \item{sp.effects}{if \code{sp.obj}
  specifies a spatial model without pre-calculated spatial effects then \code{spDiag} calculates the
    spatial effects. }
}

\seealso{\code{\link{bayesLMRef}}, \code{\link{bayesLMConjugate}},
  \code{\link{spGGT}}, \code{\link{bayesGeostatExact}},
  \code{\link{spLM}}, \code{\link{spGLM}}, \code{\link{spMvLM}}, \code{\link{mvLM}}}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical
  modeling and analysis for spatial data. Chapman and Hall/CRC Press,
  Boca Raton,Fla.
  
  Gelfand A.E. and Ghosh, S.K. (1998). Model choice: a minimum posterior
  predictive loss approach. \emph{Biometrika}. 85:1-11.

}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}.
}

\examples{
\dontrun{
data(rf.n200.dat)

Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
w <- rf.n200.dat$w

n.samples <- 1000

##non-spatial regression
m.1 <- bayesLMRef(lm(Y~1), n.samples=n.samples)


##spatial regression
m.2 <- spLM(Y~1, coords=coords,
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=n.samples, verbose=TRUE, n.report=100)

##compare
print(spDiag(m.1))
print(spDiag(m.2, start=500, thin=2))
}
}

\keyword{misc}
