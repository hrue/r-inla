context("test 'latent besag2'")

test_that("Case 1", {
    data(Oral)
    g = system.file("demodata/germany.graph", package="INLA")
    ## use data Oral to estimate a spatial field in order to simulate a
    ## 'realistic' dataset.
    formula = Y ~ f(region, model="bym", graph=g)
    result = inla(formula, data = Oral, family = "poisson", E = E)
    x = result$summary.random$region$mean
    n = length(x)/2
    ## simulate two new datasets. 'a' is the weighting between the
    ## log.rel.risk:
    a = 2
    xx = x[1:n]+1
    x = c(a*xx, xx/a)
    E = c(Oral$E, Oral$E)
    N = 2*n
    y = rpois(N, lambda = E*exp(x))
    
    ## model='besag2' defines a model with length N = 2*graph->n, the
    ## first half is weighted with 'a' the other half is weighted with
    ## 1/a. here there is no unstructed terms.
    i = 1:N
    std = TRUE
    formula = y ~ f(i, model="besag2", graph=g, scale.model=std) -1
    r = inla(formula, family = "poisson", data = data.frame(E,y,i), E=E)
    expect_true(abs(r$summary.hyperpar[2, "mean"] - a) < 0.1)
})
