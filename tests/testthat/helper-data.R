make_gaussian_data <- function(n = 24) {
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n, sd = 0.5)

  list(
    y = y,
    x = x,
    coords = coords,
    starting = list(phi = 3, sigma.sq = 1, tau.sq = 0.2),
    tuning = list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.05),
    priors = list(
      phi.Unif = c(0.1, 30),
      sigma.sq.IG = c(2, 1),
      tau.sq.IG = c(2, 1)
    )
  )
}

make_binomial_data <- function(n = 24) {
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.2 + 0.8 * x))

  list(
    y = y,
    x = x,
    coords = coords,
    weights = rep(1, n),
    starting = list(beta = c(0, 0), phi = 3, sigma.sq = 1, w = 0),
    tuning = list(beta = c(0.1, 0.1), phi = 0.1, sigma.sq = 0.1, w = 0.1),
    priors = list(
      beta.Norm = list(rep(0, 2), diag(100, 2)),
      phi.Unif = c(0.1, 30),
      sigma.sq.IG = c(2, 1)
    )
  )
}

expect_finite_matrix <- function(x, nrow = NULL, ncol = NULL) {
  expect_true(is.matrix(x))
  expect_true(all(is.finite(x)))
  if (!is.null(nrow)) {
    expect_equal(nrow(x), nrow)
  }
  if (!is.null(ncol)) {
    expect_equal(ncol(x), ncol)
  }
}
