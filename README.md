# spBayes

[![CRAN status](https://www.r-pkg.org/badges/version/spBayes)](https://CRAN.R-project.org/package=spBayes)
[![R-CMD-check](https://github.com/finleya/spBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/finleya/spBayes/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/finleya/spBayes/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/finleya/spBayes/actions/workflows/pkgdown.yaml)

<img src="man/figures/logo.png" alt="spBayes package hex logo" align="left" width="130" style="margin: 0 1.2rem 0.6rem 0;">

`spBayes` fits Bayesian univariate and multivariate spatial and
spatio-temporal regression models for point-referenced data. The package
provides MCMC-based model fitting, posterior recovery of spatial random
effects, prediction at new locations, model diagnostics, and utilities for
building spatial covariance matrices.

<br clear="left"/>

## Installation

Install the CRAN release with:

```r
install.packages("spBayes")
```

Install the development version from GitHub with:

```r
remotes::install_github("finleya/spBayes")
```

## Basic Use

```r
library(spBayes)

set.seed(1)
n <- 50
coords <- cbind(runif(n), runif(n))
x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(n)

fit <- spLM(
  y ~ x,
  coords = coords,
  starting = list(phi = 3, sigma.sq = 1, tau.sq = 1),
  tuning = list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.1),
  priors = list(
    phi.Unif = c(0.1, 30),
    sigma.sq.IG = c(2, 1),
    tau.sq.IG = c(2, 1)
  ),
  cov.model = "exponential",
  n.samples = 100,
  verbose = FALSE
)

summary(fit$p.theta.samples)
```

## Functionality

The package provides:

- Gaussian spatial linear models with `spLM()`.
- Binomial and Poisson spatial generalized linear models with `spGLM()`.
- Multivariate spatial linear and generalized linear models.
- Spatially varying coefficient and dynamic spatio-temporal models.
- Predictive-process models for larger point-referenced data.
- Posterior recovery, prediction, diagnostics, and covariance utilities.

## Citations

If you use `spBayes`, please cite:

Finley, A. O., Banerjee, S., and Carlin, B. P. (2007). spBayes: An R
Package for Univariate and Multivariate Hierarchical Point-Referenced
Spatial Models. *Journal of Statistical Software*, 19(4), 1-24.

Finley, A. O., Banerjee, S., and Gelfand, A. E. (2015). spBayes for
Large Univariate and Multivariate Point-Referenced Spatio-Temporal Data
Models. *Journal of Statistical Software*, 63(13), 1-28.
doi:10.18637/jss.v063.i13.

Finley, A. O. and Banerjee, S. (2020). Bayesian spatially varying
coefficient models in the spBayes R package. *Environmental Modelling &
Software*, 125, 104608. doi:10.1016/j.envsoft.2019.104608.
