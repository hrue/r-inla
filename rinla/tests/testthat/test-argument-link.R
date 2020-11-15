context("test 'argument link'")

test_that("Case 1", {
    set.seed(12345)
    ## simple poisson regression
    n = 100
    x = runif(n)
    eta = 1 + x
    lambda = exp(eta)
    y = rpois(n, lambda = lambda)
    
    ## missing values:
    y[1:3] = NA
    y[(n-2):n] = NA
    
    ## link = 1 is a shortcut for rep(1, <n>) where <n> is the appropriate
    ## length. here '1' is a reference to the first 'family', ie
    ## 'family[1]'
    r = inla(y ~ 1 + x,  family = "poisson",
            data = data.frame(y, x), 
            control.predictor = list(compute = TRUE, link = 1))
    expect_true(all(abs(r$summary.fitted.values$mean - lambda) < 0.6))
    
    ## this is the formally correct way, defining 'link' only where there
    ## are missing values. entries in 'link' for which the observation is
    ## not NA, is currently ignored.
    link = rep(NA, n)
    link[which(is.na(y))] = 1
    r = inla(y ~ 1 + x,  family = "poisson",
            data = data.frame(y, x), 
            control.predictor = list(compute = TRUE, link = link))
    expect_true(all(abs(r$summary.fitted.values$mean - lambda) < 0.6))

    
    ## for more than one likelihood, use '2' to refer to the second
    ## likelihood. here we just split the data in two, and assign the
    ## second half the nbinomial distribution.
    n2 = n %/% 2L
    Y = matrix(NA, n, 2)
    Y[1:n2, 1] = y[1:n2]
    Y[1:n2 + n2, 2] = y[1:n2 + n2]
    link = rep(NA, n)
    link[which(is.na(y[1:n2]))] = 1
    link[n2  + which(is.na(y[1:n2 + n2]))] = 2
    
    r = inla(Y ~ 1 + x,  family = c("poisson", "nbinomial"), 
            data = list(Y=Y, x=x), 
            control.predictor = list(compute = TRUE, link = link))
    expect_true(all(abs(r$summary.fitted.values$mean - lambda) < 0.8))
})
