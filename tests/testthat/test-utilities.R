test_that("distance and covariance utilities return expected shapes", {
  coords.1 <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
  coords.2 <- matrix(c(0, 1, 1, 1, 2, 2), ncol = 2, byrow = TRUE)

  d <- iDist(coords.1, coords.2)
  expect_finite_matrix(d, nrow = 2L, ncol = 3L)
  expect_equal(d[1, 1], 1)

  cov <- mkSpCov(coords.1, K = matrix(1), Psi = matrix(0.1), theta = 3, "exponential")
  expect_finite_matrix(cov, nrow = 2L, ncol = 2L)
  expect_equal(diag(cov), rep(1.1, 2))
})

test_that("pointsInPoly returns point membership indices", {
  poly <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
  pts <- matrix(c(0.5, 0.5, 1.5, 0.5), ncol = 2, byrow = TRUE)

  inside <- pointsInPoly(poly, pts)

  expect_equal(inside, 1L)
})
