context("test 'rw1o1diid'")

test_that("rw1o1diid matches bym2 with linear chain graph", {
  set.seed(42)
  n <- 50
  time <- 1:n
  y <- cumsum(rnorm(n, sd = 0.3)) + rnorm(n, sd = 1)

  r1 <- inla(
    y ~ f(time, model = inla.rw1o1diid(n)) - 1,
    data = data.frame(y, time)
  )

  i_idx <- c(seq_len(n - 1L), seq.int(2L, n))
  j_idx <- c(seq.int(2L, n), seq_len(n - 1L))
  graph <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(n, n))
  r2 <- inla(
    y ~ f(time, model = "bym2", graph = graph, scale.model = TRUE, n = n) - 1,
    data = data.frame(y, time)
  )

  expect_equal(
    r1$summary.random[[1]]$mean,
    r2$summary.random[[1]]$mean,
    tolerance = 1e-6
  )
  expect_equal(
    r1$summary.hyperpar[, "mean"],
    r2$summary.hyperpar[, "mean"],
    tolerance = 1e-6
  )
})

test_that("rw1o1diid input validation", {
  expect_error(inla.rw1o1diid(2))
  expect_error(inla.rw1o1diid(-1))
  expect_error(inla.rw1o1diid(c(10, 20)))
})
