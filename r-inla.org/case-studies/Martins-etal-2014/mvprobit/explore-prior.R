prior.rho.exch = function(rho, p=4L, lambda=1, log=FALSE)
{
    ## > d := (rho,p) -> sqrt(-log( (1+(p-1)*rho) * (1-rho)^(p-1) ));
    ## > a:= 1/2 * lambda*exp(-lambda*d(rho,p))*abs(diff(d(rho,p), rho));
    ## > with(CodeGeneration);
    ## > C(a, 'optimize');
    
    invalid = (rho <= -1/(p-1.0) | rho >= 1.0)
    valid = (!invalid)

    rho = rho[valid]
    pow = function(a, b) a^b
    t1 = p - 1;
    t3 = rho * t1 + 1;
    t4 = 1 - rho;
    t5 = pow(t4,  t1);
    t6 = t5 * t3;
    t7 = log( t6);
    t8 = sqrt(-t7);
    t10 = exp(-t8 * lambda);
    t23 = abs(0.1e1 /  t5 /  t3 *  (t5 * t1 - 1 / t4 * t1 * t6) / t8);
    t25 = t23 * t10 * lambda / 0.4e1;

    prior = numeric(length(valid))
    prior[invalid] = (if (log) -Inf else 0.0)
    prior[valid] = (if (log) log(t25) else t25)

    return (prior)
}

explore = function(n, nsim, lambda)
{
    r= INLA:::inla.mclapply(
            1:nsim,
            function(k) {
                theta = inla.pc.cormat.rtheta(1, p=n, lambda = lambda)
                R = inla.pc.cormat.permute(inla.pc.cormat.theta2R(theta))
                return(R[1, 2])
            })
    return (unlist(r))
}

do.explore = function()
{
    n = 4
    nsim = 100000

    lambdas = seq(0.01, 1.1, by=0.05)
    k = 1:length(lambdas)
    xx = INLA:::inla.mclapply(k,
            function(k) {
                r = explore(n, nsim, lambdas[k])
                return(r)
            })

    prob = unlist(lapply(xx, function(x) sum(abs(x) > 0.5)/length(x)))
    smooth = lowess(lambdas, log(prob), f=1/10)
    cdf = splinefun(smooth$x, exp(smooth$y))
    inv.cdf = splinefun(exp(smooth$y), smooth$x)

    ## these are the values of lambda for the saturated prior model
    lambda.low = inv.cdf(0.8)
    lambda.high = inv.cdf(0.2)

    x.low = explore(n, nsim*10, lambda.low)
    x.high = explore(n, nsim*10, lambda.high)
    hist(x.low, n=100,
         xlab = expression(rho),
         ylab = "Density",
         main = "",
         prob=TRUE)
    dev.print(postscript, file="prior-low.ps")
    system("which psfix && psfix prior-low.ps")

    hist(x.high, n=100,
         xlab = expression(rho),
         ylab = "Density",
         main = "",
         prob=TRUE)
    dev.print(postscript, file="prior-high.ps")
    system("which psfix && psfix prior-high.ps")

    ## and we plot the exchangable prior for the same lambda-values
    rho = 2/(1+exp(-seq(-20, 20, len=10000)))-1
    plot(rho, prior.rho.exch(rho, p=n, lambda = lambda.low), type="l", lwd=2,
         bty = "l", 
         xlim = c(-1/(n-1), 1),
         ylim = c(0, 5),
         xlab = expression(rho),
         ylab = "Density")
    lines(rho, prior.rho.exch(rho, p=n, lambda = lambda.high), lty=2, lwd=2)
    dev.print(postscript, file="mvprobit-prior-exch.ps")
    system("which psfix && psfix mvprobit-prior-exch.ps")
    
    ## this gives lambda.low=0.1 and lambda.high=1.0 approximately.
    ## [1] "lambda.low 0.0990822812004724"
    ## [1] "lambda.high 0.998264575560278"

    print(paste("lambda.low",  lambda.low))
    print(paste("lambda.high",  lambda.high))
}
