context("test 'likelihood cbinomial'")

test_that("Case 1", {
    set.seed(123)
    N=10000
    a = -1
    b = 1
    z = rnorm(N, sd=0.1)
    eta = a + b*z
    n = sample(c(1, 2, 3, 4), size=N, replace=TRUE)
    p = exp(eta)/(1 + exp(eta))
    prob = 1.0 - (1-p)^n
    k = sample(c(1,5,10,15), size=N, replace=TRUE)
    y = rbinom(N, size=k, prob = prob)
    data = list(y=y,z=z)
    formula = y ~ 1+z
    result = inla(formula, family = "cbinomial", data = data,
            Ntrials=cbind(k, n))
    expect_true(abs(result$summary.fixed[1,"mean"] -a) < 0.05)
    expect_true(abs(result$summary.fixed[2,"mean"] -b) < 0.05)
})   
