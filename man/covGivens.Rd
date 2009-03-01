\name{covGivens}
\alias{covGivens}          
\title{Returns covariance matrix C given components of P and L in C=PLP'}

\description{
  Given components of \eqn{P}{P} and \eqn{\Lambda}{L} in \eqn{C=P\Lambda
    P^T}{C=PLP'}, \code{covGivens} returns \eqn{C}{C}.
}

\usage{
 covGivens(theta, lambda, p, strict=TRUE)
}

\arguments{
  \item{theta}{a p(p-1)/2 vector of Givens angles that are
    restricted to lie in the interval (\eqn{-\pi/2,\pi/2}{-pi/2,pi/2}) for uniqueness and
    positive definiteness of the covariance matrix. These are arranged
    in column major order in P. }
  \item{lambda}{a vector of length p of positive eigenvalues which form
  the diagonal matrix \eqn{\Lambda}{L}.}
  \item{p}{refers to the dimension of the desired \eqn{p\times p}{pxp}
    covariance matrix.  }
  \item{strict}{checks that theta elements lie in the interval
    (\eqn{-\pi/2,\pi/2}{-pi/2,pi/2}) and the lambda elements are
  positive. }
}

\value{
  A list with the following tags:
  \item{C}{\eqn{p\times p}{pxp} covariance matrix.}
  \item{P}{\eqn{p\times p}{pxp} P matrix.}
  \item{Lambda}{diagonal elements of the \eqn{p\times p}{pxp} \eqn{\Lambda}{L}
  matrix. }
}

\details{
  The spectral decomposition of a matrix C is given by
  \eqn{P\Lambda P^T}{PLP'}, where P is an orthogonal matrix of
  eigenvectors and \eqn{\Lambda}{L} is a diagonal matrix of positive eigenvalues. In the case of distinct
  eigenvalues, \eqn{P}{P} can be expressed as \eqn{P = G_{12}, G_{13},
    G_{23}, ..., G_{p-1,p}}{P = G_{12}, G_{13}, G_{23}, ..., G_{p-1,p}}
  where \eqn{G_{ij}}{G_{i,j}} is a
  \eqn{p\times p}{pxp} matrix with \eqn{cos(\theta_{tl})}{cos(\theta_{tl})} in the \emph{t}th and \emph{l}th diagonal elements,
  \eqn{sin(\theta_{tl})}{sin(\theta_{tl})} in the \emph{tl}th and \eqn{-sin(\theta_{tl})}{-sin(\theta_{tl})} in the \emph{lt}th elements, zeros on the rest of the off-diagonal elements and 1's on the rest of
  the diagonal. The p(p-1)/2 angles are called Givens angles and are restricted to lie in the interval (\eqn{-\pi/2,\pi/2}{-pi/2,pi/2}) for uniqueness and
  positive definiteness of the resulting covariance matrix.
}

\references{
  Daniels, M.J and R.E. Kass (1999) Nonconjugate Bayesian estimation of covariance matrices and its use in hierarchical
  models, J. Amer. Statist. Assoc. 94:1254-1263.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
}

\examples{

\dontrun{
########################################
##Syntetic covariance matrix IIIR_1 from
##Daniels, M.J. and R.E. Kass (1999)
##JASA 94(448):1254-1263.
########################################
p <- 5
theta <- rep(pi/4, p*(p-1)/2)
lambda <- c(1, 0.75, 0.56, 0.06, 0.003)

cov <- covGivens(theta, lambda, p)$C
cov2cor(cov)
}
}

\keyword{misc}
