context("test 'likelihood poisson'")

test_that("Case 1", {
    set.seed(123)
    y <- 4
    Expect <- 0.25
    a <- 0; stopifnot(a==0)
    b <- 1.38
    pdat <- list(y = y, E = Expect)
    mod.pois <- inla(y ~ 1, 
                     family='poisson',
                     E=E,
                     data=pdat,
                     control.compute=list(mlik=T),
                     control.fixed=list(mean.intercept=a, prec.intercept=1/b^2),
                     control.inla = list(strategy = "laplace"))
    m <- mod.pois$marginals.fixed[[1]]
    t12 <- inla.emarginal(function(x) c(x, x^2), m)
    t1 <- t12[1]
    v1 <- t12[2]-t12[1]^2
    fac = 10
    xx = seq(t1 - fac*sqrt(v1), t1 + fac*sqrt(v1), length = 10000)
    dx = diff(xx)[1]
    ## Note: working on log-scale to avoid potential numerical issues
    ff = dpois(y, lambda = Expect * exp(xx), log = TRUE) +
        dnorm(xx, mean = a,  sd = b, log = TRUE)
    ff = ff - max(ff)
    ff = exp(ff) / sum(exp(ff)) / dx
    mm = inla.smarginal(list(x = xx, y = ff))
    tt12 <- inla.emarginal(function(x) c(x, x^2), mm)
    tt1 <- tt12[1]
    vv1 <- tt12[2]-tt12[1]^2

    expect_equal(t1, tt1, tolerance = 1e-3)
    expect_equal(v1, vv1, tolerance = 1e-1)
    
    integral = integrate(
            function(x, y, Expect, a, b)
            dnorm(x, mean=a, sd = b) * dpois(y, lambda = Expect*exp(x)),
            lower = t1 - fac*sqrt(v1), upper = t1 + fac*sqrt(v1),
            y=y, Expect=Expect, a=a, b=b)
    true.mlik = log(integral$value)

    expect_equal(mod.pois$mlik[[1]], true.mlik, tolerance = 1e-2)
})
