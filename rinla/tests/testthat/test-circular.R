context("test 'likelihood circular'")

test_that("Case 1", {
    set.seed(123)
    n = 10000
    z = rnorm(n, sd=0.3)
    eta = 1 + z
    y.pred = inla.link.invtan(eta)
    
    kappa = 10
    x = seq(-pi, pi, len = 10000)
    d = exp(kappa*cos(x))
    dd = cumsum(d)
    dd = dd /max(dd)
    cn.icdf.func = splinefun(dd, x, method = "monoH.FC")
    rcn = function(n) cn.icdf.func(runif(n))
    y = y.pred + rcn(n)
    
    formula = y ~ 1 + z
    r = inla(formula,  data = data.frame(y, z),
            family = "circularnormal",
            control.inla = list(cmin = -Inf))
    expect_true(all(abs(r$summary.fixed[, "mean"] - c(1, 1)) < 0.05))
    expect_true(all(abs(r$summary.hyperpar[1, "mean"] - kappa) < 0.1))
})
