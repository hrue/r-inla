context("test 'likelihood beta'")

test_that("Case 1", {
    set.seed(123)
    ## the precision parameter in the beta distribution
    phi = 5
    
    ## generate simulated data
    n = 1000
    z = rnorm(n, sd=0.2)
    eta = 1 + z
    mu = exp(eta)/(1+exp(eta))
    a = mu * phi
    b = -mu * phi + phi
    y = rbeta(n, a, b)
    
    ## estimate the model
    formula = y ~ 1 + z
    r = inla(formula, data = data.frame(y, z), family = "beta",
            control.compute = list(cpo=TRUE))

    expect_true(abs(r$summary.hyperpar[1, "mean"] - phi) < 0.5)
})


