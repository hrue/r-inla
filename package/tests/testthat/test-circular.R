context("test 'likelihood circular'")

test_that("Case 1", {
    inla_env <- inla.get.inlaEnv()
    testthat_circumnormal <- inla_env[["enable.model.likelihood.circularnormal"]]
    withr::defer(
        inla_env[["enable.model.likelihood.circularnormal"]] <- testthat_circumnormal
    )
    inla_env[["enable.model.likelihood.circularnormal"]] <- TRUE

    set.seed(123)
    n = 10000
    z = rnorm(n, sd=0.3)
    eta = 1 + z
    y.pred = inla.link.invtan(eta)
    
    kappa = 10
    x = seq(-pi, pi, length.out = 10000)
    d = exp(kappa*cos(x))
    dd = cumsum(d)
    dd = dd /max(dd)
    cn.icdf.func = splinefun(dd, x, method = "monoH.FC")
    rcn = function(n) cn.icdf.func(runif(n))
    y = y.pred + rcn(n)
    
    formula = y ~ 1 + z
    ## 2025-07-21: There is a 0 == 1 assertion that hardcodes the disabling
    ## of this model. If it ever gets reenabled, reactivate these tests:
    expect_error(
        {r = inla(formula,  data = data.frame(y, z),
                 family = "circularnormal",
                 control.inla = list(cmin = -Inf))},
        "The inla-program exited with an error",
        fixed = TRUE
    )
    ## expect_true(all(abs(r$summary.fixed[, "mean"] - c(1, 1)) < 0.05))
    ## expect_true(all(abs(r$summary.hyperpar[1, "mean"] - kappa) < 0.1))
})
