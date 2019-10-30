## In this example we do linear regression using 'rgeneric'.
## The regression model is y = a + b*x + noise,  and we
## define 'a + b*x + tiny.noise' as a latent model.
## The dimension is length(x) and number of hyperparameters
## is 2 ('a' and 'b').
## This is a prototype example how similar sitations 
## could be approached,  where essentially the latent model is a
## model for the 'mean' only.

rgeneric.linear.regression =
    function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", 
                     "log.prior", "quit"),
             theta = NULL)
{  
    ## the environment of this function
    envir = parent.env(environment())
  
    ## artifical high precision to be added to the mean-model
    prec.high = exp(15)
    
    interpret.theta = function() {
        return(list(a = theta[1L], b = theta[2L]))
    }
    
    graph = function() {
        G = Diagonal(n = length(x), x=1)
        return(G)
    } 
    
    Q = function() {
        Q = prec.high * graph()
        return(Q)
    }
    
    mu = function() {
        par = interpret.theta()
        return(par$a + par$b * x)
    }

    log.norm.const = function() {
        ## the easiest is to let INLA compute this
        return(numeric(0))
    }

    log.prior = function() {
        par = interpret.theta()
        val = (dnorm(par$a, mean=0, sd = sqrt(1/0.001), log=TRUE) +
               dnorm(par$b, mean = 0, sd = sqrt(1/0.001), log=TRUE))
        return(val)
    }

    initial = function() {
        return(rep(0, 2))
    }
   
    quit = function() {
        return(invisible())
    }

    val = do.call(match.arg(cmd), args = list())
    return(val)
}

a = 1
b = 2
n = 100
x = rnorm(n)
eta = a + b*x
y = eta + rnorm(n)

rgen = inla.rgeneric.define(model = rgeneric.linear.regression, x = x)
r = inla(y ~ -1 + f(idx, model=rgen),
         data = data.frame(y, idx = 1:n))
rr = inla(y ~ 1 + x,
          data = data.frame(y, x),
          control.fixed = list(prec.intercept = 0.001, prec = 0.001))

## compare the results with the 'truth'
par(mfrow=c(2, 1))
plot(r$marginals.hyperpar[['Theta1 for idx']], type="l",  lwd=2, col="red", 
     main = "Posterior for the intercept (red=rgeneric, blue=default)")
lines(rr$marginals.fixed$'(Intercept)', lwd=2, col="blue")

plot(r$marginals.hyperpar[['Theta2 for idx']], type="l",  lwd=2, col="red", 
     main = "Posterior for the slope (red=rgeneric,  blue=default)")
lines(rr$marginals.fixed$'x', lwd=2, col="blue")
