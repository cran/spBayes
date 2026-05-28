test_that("spGLM fits a small binomial spatial model", {
  set.seed(4)
  dat <- make_binomial_data(n = 24)

  capture.output(
    fit <- spGLM(
      dat$y ~ dat$x,
      family = "binomial",
      coords = dat$coords,
      weights = dat$weights,
      starting = dat$starting,
      tuning = dat$tuning,
      priors = dat$priors,
      amcmc = list(n.batch = 3, batch.length = 4, accept.rate = 0.43),
      cov.model = "exponential",
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "spGLM")
  expect_equal(fit$family, "binomial")
  expect_equal(fit$cov.model, "exponential")
  expect_equal(dim(fit$p.beta.theta.samples), c(12L, 4L))
  expect_finite_matrix(fit$p.w.samples, nrow = 24L, ncol = 12L)
  expect_equal(dim(fit$acceptance), c(4L, 3L))
})

test_that("spGLM fits a nonspatial binomial model without coords", {
  set.seed(5)
  dat <- make_binomial_data(n = 24)

  capture.output(
    fit <- spGLM(
      dat$y ~ dat$x,
      family = "binomial",
      weights = dat$weights,
      starting = list(beta = c(0, 0)),
      tuning = list(beta = c(0.1, 0.1)),
      priors = list(beta.Norm = list(rep(0, 2), diag(100, 2))),
      amcmc = list(n.batch = 3, batch.length = 4, accept.rate = 0.43),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "nonSpGLM")
  expect_equal(fit$family, "binomial")
  expect_equal(dim(fit$p.beta.samples), c(12L, 2L))
})

test_that("spGLM validates supported families", {
  set.seed(6)
  dat <- make_binomial_data(n = 20)

  expect_error(
    spGLM(
      dat$y ~ dat$x,
      family = "gaussian",
      weights = dat$weights,
      n.samples = 5,
      verbose = FALSE
    ),
    "family must be binomial or poisson"
  )
})
