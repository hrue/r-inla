library(gsl)

inla.pc.ar.lambda = function(a = 0.5, b = 0.8, p = 4)
{
    inla.pc.ar.solve.lambda = function(pred.err.factors, nseq = 1000L) {
        ## pred.err.factor = E(1-rho^2)
        pred.err = function(lambda) {
            return(0.5 * lambda * sqrt(pi) * exp(lambda^2/4 + log_erfc(lambda/2)))
        }

        ## find lower and upper limit of lambda for the pred.err.factors given
        lambda.min = lambda.max = 1
        val.max = max(pred.err.factors)
        val.min = min(pred.err.factors)
        while(pred.err(lambda.min) > val.min) {
            lambda.min = lambda.min / 2.0
        }
        while(pred.err(lambda.max) < val.max) {
            lambda.max = lambda.max * 2.0
        }

        lambda = lambda.min * exp(seq(0, log(lambda.max/lambda.min), length = nseq))
        fun = splinefun(pred.err(lambda), lambda,  method = "monoH.FC")
        lambdas = fun(pred.err.factors)

        return (lambdas)
    }

    pred.err.factors = 1.0 - (1.0-a)*b^(0:(p-1))
    lambda = inla.pc.ar.solve.lambda(pred.err.factors)

    return (lambda)
}

inla.pc.ar.test1 = function(p=1, p.est = p+1, n = 100)
{
    pacf = c(runif(p), rep(0, p.est-p))
    phi = inla.ar.pacf2phi(pacf)
    sd.x = sqrt(1/prod(1-pacf^2))
    x = arima.sim(n = n, model = list(ar = phi))/sd.x
    lambda = c(inla.pc.ar.lambda(p = p.est, b = 0.5), rep(1, 10))
    initial = c(inla.models()$latent$ar$hyper$theta2$to.theta(pacf), rep(0, 10))
    r = (inla(
        x ~ -1 +
            f(time,  model = "ar",  order = p.est, 
              hyper = list(
                  prec = list(param = c(3, 0.01), initial = 0), 
                  pacf1 = list(param = c(lambda[1], 0), initial = initial[1]),
                  pacf2 = list(param = c(lambda[2], 0), initial = initial[2]),
                  pacf3 = list(param = c(lambda[3], 0), initial = initial[3]),
                  pacf4 = list(param = c(lambda[4], 0), initial = initial[4]),
                  pacf5 = list(param = c(lambda[5], 0), initial = initial[5]),
                  pacf6 = list(param = c(lambda[6], 0), initial = initial[6]),
                  pacf7 = list(param = c(lambda[7], 0), initial = initial[7]),
                  pacf8 = list(param = c(lambda[8], 0), initial = initial[8]),
                  pacf9 = list(param = c(lambda[9], 0), initial = initial[9]),
                  pacf10 = list(param = c(lambda[10], 0), initial = initial[10]))), 
        data = data.frame(x=x, time = 1:length(x)),
        control.family = list(hyper = list(prec = list(initial = 12,
                                                       fixed=TRUE)))))
    result = cbind(est = r$summary.hyperpar$mean[-1], true = pacf)
    inside = c()
    for(i in 1:p.est) {
        int = inla.hpdmarginal(0.95, r$marginals.hyperpar[[i+1]])
        inside[i] = (int[1] < pacf[i] && pacf[i] < int[2])
    }
    result = cbind(result, coverage = inside)
    print(round(result, digits=3))
    
    result = (cbind(acf.est = inla.ar.pacf2acf(r$summary.hyperpar$mean[-1], lag.max = 10),
                    acf.emp = c(acf(x, lag.max=10, plot=FALSE)$acf), 
                    acf.true = inla.ar.pacf2acf(pacf, lag.max=10)))
    print(round(result, digits=3))

    return (invisible())
}

