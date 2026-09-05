context("test 'likelihood poissonlognormal'")

## Case 1: intercept recovery
## Simulate y_i ~ Poisson(exp(beta0 + u_i)), u_i ~ N(0, 1/prec_true).
## Check that the posterior mean of the intercept is close to the truth.
test_that("Case 1: intercept recovery", {
    set.seed(123)
    n         <- 500
    beta0     <- 1.0
    sigma     <- 1.0
    y         <- rpois(n, exp(beta0 + rnorm(n, 0, sigma)))
    r         <- inla(y ~ 1,
                      family          = "poissonlognormal",
                      data            = data.frame(y = y),
                      control.compute = list(dic = TRUE))
    expect_true(abs(r$summary.fixed[1L, "mean"] - beta0) < 0.05)
})

## Case 2: precision recovery
## Same generative model; check that the posterior mean of the
## lognormal precision hyperparameter is close to 1/sigma^2.
test_that("Case 2: precision recovery", {
    set.seed(123)
    n         <- 500
    beta0     <- 1.0
    sigma     <- 1.0
    prec_true <- 1 / sigma^2
    y         <- rpois(n, exp(beta0 + rnorm(n, 0, sigma)))
    r         <- inla(y ~ 1,
                      family          = "poissonlognormal",
                      data            = data.frame(y = y),
                      control.compute = list(dic = TRUE))
    expect_true(abs(r$summary.hyperpar[1L, "mean"] - prec_true) < 0.3)
})

## Case 3: marginal likelihood
## For n=1 with fixed precision, INLA's log-marginal-likelihood must match
## the true value obtained by numerically integrating out the intercept.
##
## True p(y) = integral_eta  p_PLN(y | eta, prec_fixed) * N(eta; a, b^2) d_eta
##
## where p_PLN(y | eta, prec_fixed)
##       = integral_u  Poisson(y; exp(eta+u)) * N(u; 0, 1/prec_fixed) du
test_that("Case 3: marginal likelihood vs numerical integral", {
    y_obs      <- 3L
    E_obs      <- 1.0
    prec_fixed <- 1.0
    a          <- 0.0          ## prior mean on intercept
    b          <- 1.38         ## prior sd on intercept (matches test-poisson.R)

    ## true marginal likelihood by double numerical integration
    sd_fixed <- 1 / sqrt(prec_fixed)
    pln_lik  <- function(eta) {
        vapply(eta, function(e) {
            integrate(
                function(u) dpois(y_obs, E_obs * exp(e + u)) * dnorm(u, 0, sd_fixed),
                lower = -20, upper = 20
            )$value
        }, numeric(1L))
    }
    fac       <- 15
    eta_lo    <- a - fac * b
    eta_hi    <- a + fac * b
    true_mlik <- log(
        integrate(function(eta) pln_lik(eta) * dnorm(eta, a, b),
                  lower = eta_lo, upper = eta_hi)$value
    )

    ## INLA marginal likelihood with the same fixed precision and prior
    r <- inla(y ~ 1,
              family          = "poissonlognormal",
              data            = data.frame(y = y_obs),
              control.compute = list(mlik = TRUE),
              control.fixed   = list(mean.intercept  = a,
                                     prec.intercept  = 1 / b^2),
              control.family  = list(hyper = list(
                  prec = list(fixed = TRUE, initial = log(prec_fixed))
              )),
              control.inla    = list(strategy = "laplace"))

    expect_equal(r$mlik[[1L]], true_mlik, tolerance = 1e-2)
})
