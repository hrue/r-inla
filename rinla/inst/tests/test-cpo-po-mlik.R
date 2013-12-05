context("test 'cpo-po-mlik'")

test_that("Case 1", {
    set.seed(12345)
    library(INLA)
    library(mvtnorm)
    N = rpois(1, lambda = 100)
    rho = runif(1)
    sigma.eps = runif(1)
    tau.eps=1/sigma.eps^2
    sigma.x = runif(1) ## marginal stdev!
    tau.x=1/sigma.x^2  ## marginal prec!
    C0 = matrix(0, nrow=N, ncol=N)
    for(i in 1:(N-1))
    {
        C0[i,i+1]=-rho
        C0[i+1,i]=-rho
    }
    diag(C0) = 1+rho^2
    C0[1,1]=1
    C0[N,N]=1
    ## Cov for y
    Sigma = solve(tau.x * C0/(1-rho^2)) + 1/tau.eps * diag(N)

    x = rep(NA,N)
    x[1] = rnorm(1, mean = 0, sd = sigma.x*sqrt(1-rho^2))
    for(i in 2:N) {
        x[i] = rho*x[i-1] + sigma.x*rnorm(1)
    }
    y = x + sigma.eps*rnorm(N)

    d = data.frame(id=1:N, y=y)
    fit = inla(y ~ -1 +
            f(id,
              model="generic0",
              Cmatrix=C0/(1-rho^2),
              hyper=list(prec=list(initial=log(tau.x),fixed=TRUE))),
            data=d,
            control.family=list(
                    hyper=list(
                            prec=list(initial=log(tau.eps),fixed=TRUE))),
            control.compute=list(cpo=TRUE, po=TRUE))

    theta1=log(tau.x)
    theta2=log((1+rho)/(1-rho))
    fit2 = inla(y~ -1 +
            f(id,
              model="ar1",
              hyper=list(
                      theta1=list(initial=theta1,fixed=TRUE),
                      theta2=list(initial=theta2,fixed=TRUE))),
            data=d,
            control.family=list(
                    hyper=list(
                            prec=list(initial=log(tau.eps),fixed=TRUE))),
            control.compute=list(cpo=TRUE, po=TRUE))

    ##Checking marginal likelihood, seems ok
    dens = dmvnorm(y,rep(0,N),Sigma,log=TRUE)
    err.mlik = sum(abs(c(fit$mlik[1]+
            0.5*determinant(C0/(1-rho^2),log=TRUE)$mod,
            fit2$mlik[2]) -  dens))

    ##Checking cpo
    cpo = numeric(N)
    for(i in 1:N) {
        D = diag(N)
        D[i, i] = 0
        QQ = tau.x * C0/(1-rho^2) + tau.eps * D
        bb = tau.eps * y
        bb[i] = 0
        sQQ = solve(QQ)
        mu = (sQQ %*% bb)[i]
        var = sQQ[i, i] + 1/tau.eps
        cpo[i] = dnorm(y[i], mean = mu, sd = sqrt(var))
    }

    ##po
    D = diag(N)
    QQ = tau.x * C0/(1-rho^2) + tau.eps * D
    bb = tau.eps * y
    sQQ = solve(QQ)
    mu = c(sQQ %*% bb)
    var = diag(sQQ) + 1/tau.eps
    po = dnorm(y, mean = mu, sd = sqrt(var))

    res.po = cbind(fit = fit$po$po, fit2 = fit2$po$po)
    err.po = mean(abs(res.po - po))
    res = cbind(fit = fit$cpo$cpo, fit2 = fit2$cpo$cpo)
    err.cpo = mean(abs(res -cpo))

    expect_true(err.mlik < 0.001)
    expect_true(err.po < 0.00001)
    expect_true(err.cpo < 0.00001)
})
