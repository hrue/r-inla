context("test 'likelihood betabinomial'")

test_that("Case 1", {
    set.seed(123)
    ## overdispersion parameter in the betabinomial
    rho = 0.1
    pzero = 0.2
    type=0
    n = 1000
    z = rnorm(n, sd=0.2)
    Ntrials = sample(1:20, n, replace=TRUE)
    eta = 1 + z
    p.eta = exp(eta)/(1+exp(eta))
    a = p.eta * (1-rho)/rho
    b = (p.eta * rho - p.eta - rho + 1)/rho
    p = rbeta(n, a, b)
    
    y = numeric(n)
    for(i in 1:n) {
        pp = dbinom(0:Ntrials[i], size = Ntrials[i], prob = p[i])
        if (type == 0) {
            pp[-1] = (1-pzero) * pp[-1]/(1-pp[1])
            pp[1] = pzero
        } else {
            pp = pp * (1-pzero)
            pp[1] = pp[1] + pzero
        }
        y[i] = sample(0:Ntrials[i], 1, prob = pp)
    }
    formula = y ~ 1 + z
    data = data.frame(y, z)
    r = inla(formula, data = data,
            family = paste("zeroinflatedbetabinomial", type, sep=""),
            Ntrials=Ntrials, 
            control.compute=list(cpo=TRUE))

    expect_true(abs(r$summary.hyperpar[1, "mean"] - rho) < 0.05)
    expect_true(abs(r$summary.hyperpar[2, "mean"] - pzero) < 0.025)
})

