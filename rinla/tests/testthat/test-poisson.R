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
    ff = dpois(y, lambda = Expect * exp(xx)) * dnorm(xx, mean = a,  sd = b)
    ff = ff/max(ff)
    ff = ff / sum(ff) / dx
    mm = inla.smarginal(list(x = xx, y = ff))
    tt12 <- inla.emarginal(function(x) c(x, x^2), mm)
    tt1 <- tt12[1]
    vv1 <- tt12[2]-tt12[1]^2

    expect_true(abs(t1 - tt1) < 0.01)
    expect_true(abs(v1 - vv1) < 0.01)
    
    integral = integrate(
            function(x, y, Expect, a, b)
            dnorm(x, mean=a, sd = b) * dpois(y, lambda = Expect*exp(x)),
            lower = t1 - fac*sqrt(v1), upper = t1 + fac*sqrt(v1),
            y=y, Expect=Expect, a=a, b=b)
    true.mlik = log(integral$value)

    expect_true(abs(true.mlik - mod.pois$mlik[[1]]) < 0.02)
})
