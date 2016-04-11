## translated and modified C-code from Leo, leo-code.c

aux.pow = function(a, b) return(a^b)
aux.sqr = function(a) return(a * a)

aux.trunc.norm.left = function(a)
{
    if (a < 0) {
        stop = 0
        while (stop == 0) {
            x = rnorm(1)
            stop = (x > a)
        }
    } else {
        a_star = 0.5 * (a + sqrt(a * a + 4))
        stop = 0
        while (stop == 0) {
            x = a - log(runif(1))/a_star
            stop = ((log(runif(1)) < (x * (a_star - 0.5 * x) - aux.sqr(a_star) * 
                0.5)))
        }
    }
    return(x)
}

## generates a sample from the logistic distribution with mean mean truncated to
## values only left from zero

aux.trunc.logistic.left = function(mean)
{
    Fa = inla.link.invlogit(mean)
    arg = Fa + runif(1) * (1 - Fa)
    result = mean - inla.link.logit(arg)
    return(result)
}

## left = 1: sampling left from zero left = 0: sampling right from zero
aux.trunc.logistic = function(mean, left)
{
    if (left == 1) {
        return(aux.trunc.logistic.left(mean))
    } else if (left == 0) {
        return(-aux.trunc.logistic.left(-mean))
    } else {
        stop()
    }
    return(0)
}

aux.update.z = function(eta, y)
{
    n = length(eta)
    z = numeric(n)
    for (i in 1:n) {
        z[i] = aux.trunc.logistic(eta[i], 1 - y[i])
    }
    return(z)
}

aux.update.aux.variance = function(z, eta)
{
    n = length(z)
    lambda = numeric(n)
    for (i in 1:n) {
        lambda[i] = aux.lambda.fc((z[i] - eta[i])^2)
    }
    
    return(lambda)
}

## Generator for inverse Gaussian distribution taken from Devroye, page 149
aux.IG = function(mu, lambda)
{
    mu2 = aux.sqr(mu)
    N = rnorm(1)
    Y = aux.sqr(N)
    X1 = mu + mu2 * Y/(2 * lambda) - mu/(2 * lambda) * sqrt(4 * mu * lambda * Y + 
        mu2 * aux.sqr(Y))
    if (runif(1) <= mu/(mu + X1)) {
        return(X1)
    } else {
        return(mu2/X1)
    }
}

aux.GIG = function(lambda, psi, chi)
{
    ## returns a sample from a generalized inverse gaussian distribution with
    ## parameters lambda, chi, psi Devroye, page 478 parametrization assume lambda > 0
    
    if ((lambda == 0.5) && (psi == 1)) {
        ## VERSION 1: LEO 14 arithmetic operations return r/aux.IG(1., r)
        
        r = sqrt(chi)
        N = rnorm(1)
        Y = aux.sqr(N)
        ## take into account that r can be small. fixed by hrue
        if (r/Y < 1e-07) {
            X1 = r/Y
        } else {
            X1 = 1 + (Y - sqrt((4 * r + Y) * Y))/(2 * r)
        }
        if (runif(1) <= 1/(1 + X1)) {
            return(r/X1)
        } else {
            return(r * X1)
        }
    } else {
        ## this first part of the program uses a ratio-of-unforms method to draw a sample
        ## from the generalised-inverse-gaussian density.
        
        quart = 0.25
        one = 1
        vsmall = 1e-29
        vlarge = 1e+31
        
        h = lambda
        b = sqrt(psi * chi)
        
        if ((h < 0) || (b < 0)) 
            stop()
        if (1) {
            if (h > quart * b * sqrt(vlarge)) 
                stop()
            e = b * b
            d = h + one
            ym = (-d + sqrt(d * d + e))/b
            if (ym < vsmall) 
                stop()
            d = h - one
            xm = (d + sqrt(d * d + e))/b
            d = 0.5 * d
            e = -quart * b
            r = xm + one/xm
            w = xm * ym
            a = aux.pow(w, (-0.5 * h)) * sqrt(xm/ym) * exp(-e * (r - ym - one/ym))
            if (a < vsmall) 
                stop()
            c = -d * log(xm) - e * r
        }
        
        ok = 0
        while (!ok) {
            r1 = runif(1)
            r2 = runif(1)
            x = a * r2/r1
            if (log(r1) < (d * log(x) + e * (x + one/x) + c)) 
                ok = 1
        }
        fn_val = sqrt(chi/psi) * x
        return(fn_val)
    }
}

aux.f1old = function(x, j)
{
    ## first series approximation, guranteed to be monotone for x > 1.333
    a = aux.sqr(j + 1) * exp(-0.5 * x * (aux.sqr(j + 1) - 1))
    return(a)
}
aux.f2old = function(x, j)
{
    ## second series approximation, guranteed to be monotone for x < 1.333
    if ((j%%2) == 1) {
        ## odd
        a = (x/(pi^2)) * exp(-((aux.sqr(j) - 1) * (pi^2))/(2 * x))
    } else {
        ## even
        a = aux.sqr(j + 1) * exp(-((aux.sqr(j + 1) - 1) * (pi^2))/(2 * x))
    }
    return(a)
}

## returns sample from full conditional of lambda
aux.lambda.fc = function(chi)
{
    ok_out = 0
    while (!ok_out) {
        fn_val = aux.GIG(0.5, 1, chi)
        u = runif(1) * 4
        if (fn_val > 1.334) {
            j = 1
            factor = 4
            upper = 1
            ## apply squeezing
            while (1) {
                ## first adjust the lower bound using the alternating series f1() - given below
                lower = upper - aux.f1old(fn_val, j)
                if (u < factor * lower) {
                  ## if the draw of the acceptance prob is below the lower bound then you are
                  ## definatly a draw from the density - ACCEPT
                  ok_out = 1
                  break
                }
                ## now adjust the upper bound
                upper = lower + aux.f1old(fn_val, j + 1)
                if (u > factor * upper) {
                  ## if the draw of the acceptance prob is above the upper bound then you are
                  ## definatly NOT a draw from the density - REJECT
                  ok_out = 0
                  break
                }
                ## else u lies somewhere inbetween (lower, upper) so we're not sure and we must
                ## continue
                j = j + 2
            }
        } else {
            ## you are in the other region (the above series f1() is not guarenteed to be
            ## monotone hence we must find other monotone series
            j = 1
            
            ## aa is simply the supremum under the rejection sampler which lies at the
            ## boundary
            aa = 0.5 * log(2 * pi) + 2 * log(pi) + log(4) - 2.5 * log(fn_val) - (aux.sqr(pi)/(2 * 
                fn_val)) + 0.5 * fn_val
            factor = exp(aa)
            upper = 1
            while (1) {
                ## this bit is the same as above but we use the series f2()
                lower = upper - aux.f2old(fn_val, j)
                if (u < factor * lower) {
                  ok_out = 1
                  break
                }
                upper = lower + aux.f2old(fn_val, j + 1)
                if (u > factor * upper) {
                  ok_out = 0
                  break
                }
                j = j + 2
            }
        }
    }
    
    return(fn_val)
}
 
