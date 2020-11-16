context("test 'latent ar'")

test_that("Case 1", {
    for(i in 1:20) {
        pac = runif(i)
        phi = inla.ar.pacf2phi(pac)
        pac2 = inla.ar.phi2pacf(phi)
        expect_true( mean(abs(pac2-pac)) < 0.00001)
    }
})

test_that("Case 2", {
    set.seed(123)
    n = 300
    phi = 0.9
    y = arima.sim(n, model = list(ar = phi))
    y = scale(y)*2
    idx = 1:n
    param.phi = c(0, 1)
    param.prec = c(1, 0.01)
    r1 = inla(y ~ -1 + f(idx, model='ar1',
            hyper = list(
                    prec = list(initial = 1, fixed=FALSE, param = param.prec),
                    rho = list(initial = 3, fixed=FALSE, param = param.phi))),
            family = "gaussian", 
            control.family = list(initial = 5, fixed=TRUE), 
            data = data.frame(y, idx))

    r = inla(y ~ -1 + f(idx, model='ar',
            order = 1, 
            hyper = list(
                    prec = list(initial = 1, fixed=FALSE, param = param.prec),
                    theta2 = list(initial = 3, fixed=FALSE, param = param.phi))), 
            family = "gaussian", 
            control.family = list(initial = 5, fixed=TRUE), 
            data = data.frame(y, idx)) 
    expect_true(all(abs(r1$summary.hyperpar[, "mean"] - r$summary.hyperpar[, "mean"]) < 0.001))
})
    
test_that("Case 3", {
    set.seed(1234)
    n = 1000
    p = 3
    pacf = runif(p)
    phi = inla.ar.pacf2phi(pacf)
    y = arima.sim(n, model = list(ar = phi))
    two = 2
    y = scale(y)*two
    idx = 1:n
    param.phi = c(rep(0, p),  diag(p))
    param.prec = c(1, 0.01)
    r = inla(y ~ -1 + f(idx, model='ar', order = p), 
            family = "gaussian", 
            control.family = list(initial = 8, fixed=TRUE), 
            data = data.frame(y, idx))
    expect_true(abs(r$summary.hyperpar[1, "mode"] - 1/two^2) < 0.01)
    expect_true(all(abs(r$summary.hyperpar[-1, "mode"] - pacf) < 0.1))
})

                
                
