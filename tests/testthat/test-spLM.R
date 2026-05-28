test_that("spLM fits a small Gaussian spatial model and supports recovery/prediction", {
  set.seed(1)
  dat <- make_gaussian_data(n = 24)

  fit <- spLM(
    dat$y ~ dat$x,
    coords = dat$coords,
    starting = dat$starting,
    tuning = dat$tuning,
    priors = dat$priors,
    cov.model = "exponential",
    n.samples = 10,
    verbose = FALSE
  )

  expect_s3_class(fit, "spLM")
  expect_equal(fit$cov.model, "exponential")
  expect_equal(colnames(fit$p.theta.samples), c("sigma.sq", "tau.sq", "phi"))
  expect_equal(dim(fit$p.theta.samples), c(10L, 3L))
  expect_finite_matrix(fit$coords, nrow = 24L, ncol = 2L)

  capture.output(
    recovered <- spRecover(fit, start = 5, verbose = FALSE)
  )
  expect_s3_class(recovered, "spLM")
  expect_equal(dim(recovered$p.beta.recover.samples), c(6L, 2L))
  expect_equal(dim(recovered$p.theta.recover.samples), c(6L, 3L))
  expect_finite_matrix(recovered$p.w.recover.samples, nrow = 24L, ncol = 6L)

  capture.output(
    pred <- spPredict(
      recovered,
      pred.coords = dat$coords[1:3, ],
      pred.covars = cbind(1, dat$x[1:3]),
      start = 5,
      verbose = FALSE
    )
  )
  expect_finite_matrix(pred$p.y.predictive.samples, nrow = 3L, ncol = 6L)
})

test_that("spLM falls back to nonspatial reference regression without coords", {
  set.seed(2)
  dat <- make_gaussian_data(n = 20)

  fit <- spLM(dat$y ~ dat$x, n.samples = 8, verbose = FALSE)

  expect_s3_class(fit, "bayesLMRef")
  expect_equal(dim(fit$p.beta.tauSq.samples), c(8L, 3L))
})

test_that("spLM validates required covariance inputs", {
  set.seed(3)
  dat <- make_gaussian_data(n = 20)

  expect_error(
    spLM(
      dat$y ~ dat$x,
      coords = dat$coords,
      starting = dat$starting,
      tuning = dat$tuning,
      priors = dat$priors,
      n.samples = 5,
      verbose = FALSE
    ),
    "cov.model must be specified"
  )
})
