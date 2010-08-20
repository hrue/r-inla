

inla.sb.1 = function()
{


    ### test marginal-routines

    x = seq(-6,6, len=1000)
    y = dnorm(x)

    marginal = cbind(x,y)

    print(c(1,inla.expectation(function(x) 1, marginal)))
    print(c(0,inla.expectation(function(x) x, marginal)))
    print(c(1,inla.expectation(function(x) x^2, marginal)))
    print(c(0,inla.expectation(function(x) x^3, marginal)))
    print(c(3,inla.expectation(function(x) x^4, marginal)))

    print("qnorm")
    for (p in seq(0.0001, 0.9999, len=10)) {
        print(c(p, qnorm(p), inla.qmarginal(p, marginal)))
    }

    print("pnorm")
    for (q in seq(-5,5,len=10)) {
        print(c(q, pnorm(q), inla.pmarginal(q, marginal)))
    }
    
    print("dnorm")
    for (x in seq(-5,5,len=10)) {
        print(c(x, dnorm(x), inla.dmarginal(x, marginal)))
    }
    print("dnorm log=TRUE")
    for (x in seq(-5,5,len=10)) {
        print(c(x, dnorm(x, log=TRUE), inla.dmarginal(x, marginal, log=TRUE)))
    }

}

