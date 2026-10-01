context("test 'rw2o1diid'")

test_that("rw2o1diid constructor returns a usable rgeneric object", {
  n <- 50
  m <- inla.rw2o1diid(n)
  expect_s3_class(m, "inla.rgeneric")
  expect_equal(m$f$model, "rgeneric")
  expect_equal(m$f$n, n)

  fn <- m$f$rgeneric$definition

  init <- do.call(fn, list(cmd = "initial", theta = NULL))
  expect_length(init, 2)

  G <- do.call(fn, list(cmd = "graph", theta = NULL))
  expect_equal(dim(G), c(n, n))

  Q <- do.call(fn, list(cmd = "Q", theta = init))
  expect_equal(dim(Q), c(n, n))
  expect_true(isSymmetric(as.matrix(Q)))

  expect_true(is.finite(do.call(
    fn,
    list(cmd = "log.norm.const", theta = init)
  )))
  expect_true(is.finite(do.call(fn, list(cmd = "log.prior", theta = init))))
})

test_that("rw2o1diid Q matches the (tau, phi) mixture formula", {
  n <- 30
  m <- inla.rw2o1diid(n)
  fn <- m$f$rgeneric$definition

  R <- as.matrix(inla.rw(n, order = 2, scale.model = TRUE))

  ## phi -> 1: Q approaches tau * R
  theta_hi <- c(log(2), 20)
  Q_hi <- as.matrix(do.call(fn, list(cmd = "Q", theta = theta_hi)))
  expect_lt(max(abs(Q_hi - 2 * R)), 1e-3)

  ## phi -> 0: Q approaches tau * I
  theta_lo <- c(log(3), -20)
  Q_lo <- as.matrix(do.call(fn, list(cmd = "Q", theta = theta_lo)))
  expect_lt(max(abs(Q_lo - 3 * diag(n))), 1e-3)
})

test_that("rw2o1diid log.prior is finite, smooth, and peaks near phi = 0.5 with default prior", {
  n <- 40
  m <- inla.rw2o1diid(n)
  fn <- m$f$rgeneric$definition

  grid <- seq(-6, 6, length.out = 25)
  lps <- vapply(
    grid,
    function(t2) do.call(fn, list(cmd = "log.prior", theta = c(0, t2))),
    numeric(1)
  )
  expect_true(all(is.finite(lps)))

  ## with default (u = 0.5, alpha = 0.5), the PC prior on phi has its mode at
  ## the median phi = 0.5, i.e. logit(phi) = 0
  expect_equal(grid[which.max(lps)], 0, tolerance = 1)
})

test_that("rw2o1diid input validation", {
  expect_error(inla.rw2o1diid(4))
  expect_error(inla.rw2o1diid(-1))
  expect_error(inla.rw2o1diid(c(10, 20)))
  expect_error(inla.rw2o1diid(50, prior.tau = list(u = 1)))
  expect_error(inla.rw2o1diid(50, prior.phi = list(alpha = 0.5)))
})
