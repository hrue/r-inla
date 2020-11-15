context("test 'likelihood gamma'")

test_that("Case 1", {
    set.seed(123)
    n = 1000
    x = rnorm(n)
    eta = 1 + x
    mu = exp(eta)
    prec.scale = runif(n, min = 0.5, max = 2)
    prec.par = 1.2
    a = prec.par * prec.scale
    b = mu / (prec.par * prec.scale)
    y = rgamma(n, shape = a, scale = b)
    r = inla(y ~ 1 + x,  data = data.frame(y, x),
            scale = prec.scale, family = "gamma")

    expect_true(abs(r$summary.fixed[1, "mean"] - 1) < 0.01)
    expect_true(abs(r$summary.fixed[2, "mean"] - 1) < 0.01)
    expect_true(abs(r$summary.hyperpar[1, "mean"] - prec.par) < 0.1)
})
