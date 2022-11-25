## Export: inla.likelihood

## ! \name{inla.likelihood}
## ! \alias{likelilhood}
## !
## ! \title{Providing functions for sampling new data, evaluating pdf, cdf, and quantiles for new data.}
## !
## ! \description{This function return function to compute the pdf,cdf,quantiles, or samples for new data
## ! using the likelihood from a inla-object.}
## ! \usage{
## !     inla.likelihood(type = c(d,r,p,q),args)
## ! }
## !
## ! \arguments{
## !   \item{type}{The returned function type. The definition is similar to "rnorm","dnorm","pnorm",and "dnorm".}
## !   \item{args}{It is usually a return value from "inla.likelihood.parser", which specifies parameters, link function and transformation function of hyperparameters.}
## !}
## !\details{Details goes here}
## !\value{value goes here}
## !\author{Havard Rue \email{hrue@r-inla.org}}
## !
## !\examples{
## ! }
weibull.likelihood = function(args){
    stopifnot(args$variant == 1)
    shape = inla.models()$likelihood$weibull$hyper$theta$from.theta(args$theta)
    scale = args$inv.link.function(args$linear.predictor)
    islog = args$log
    lower.tail = args$lower.tail
    if(args$type == "r"){
        fun = function(n){
            return(rweibull(n,scale = scale,shape = shape))
        }
    }else if(args$type == "d"){
        fun = function(x){
            return(dweibull(x = x,scale = scale,shape = shape,log = islog))
        }
    }else if(args$type == "p"){
        fun = function(q){
            return(pweibull(q = q,scale = scale,shape = shape,log.p = islog,lower.tail = lower.tail))
        }
    }else if(args$type == "q"){
        fun = function(p){
            return(qweibull(p = p,scale = scale,shape = shape,log.p = islog,lower.tail = lower.tail))
        }
    }else if(args$type == "s"){
        fun = function(q){
            return(pweibull(q = q,scale = scale,shape = shape,log.p = islog,lower.tail = FALSE))
        }
    }
    
    return (fun)
}

poisson.likelihood = function(args){
    lambda = args$inv.link.function(args$linear.predictor)
    islog = args$log
    lower.tail = args$lower.tail
    if(args$type == "r"){
        fun = function(n){
            return(rpois(n = n,lambda = lambda))
        }
    }else if(args$type == "d"){
        fun = function(x){
            return(dpois(x = x,lambda = lambda,log = islog))
        }
    }else if(args$type == "p"){
        fun = function(q){
            return(ppois(q = q,lambda = lambda,log.p = islog,lower.tail = lower.tail))
        }
    }else if(args$type == "q"){
        fun = function(p){
            return(qpois(p = p,lambda = lambda,log.p  = islog,lower.tail = lower.tail))
        }
    }else if(args$type == "s"){
        fun = function(q){
            return(ppois(q = q,lambda = lambda,log.p = islog,lower.tail = FALSE))
        }
    }
    
    return (fun)
}

gaussian.likelihood = function(args){
    mean = args$inv.link.function(args$linear.predictor)
    prec = inla.models()$likelihood$gaussian$hyper$theta1$from.theta(args$theta)
    scale = args$scale
    stdev = 1/sqrt(prec*scale)
    islog = args$log
    lower.tail = args$lower.tail
    if(args$type == "r"){
        fun = function(n){
            return(rnorm(n = n,mean = mean,sd = stdev))
        }
    }else if(args$type == "d"){
        fun = function(x){
            return(dnorm(x = x,mean = mean,sd = stdev,log = islog))
        }
    }else if(args$type == "p"){
        fun = function(q){
            return(pnorm(q = q,mean = mean,sd = stdev,log.p = islog,lower.tail = lower.tail))
        }
    }else if(args$type == "q"){
        fun = function(p){
            return(qnorm(p = p,mean = mean,sd = stdev,log.p  = islog,lower.tail = lower.tail))
        }
    }else if(args$type == "s"){
        fun = function(q){
            return(pnorm(q = q,mean = mean,sd = stdev,log.p = islog,lower.tail = FALSE))
        }
    }
    
    return (fun)
}



inla.likelihood = function(type = c("d","p","r","q","s"),args){
    stopifnot(type %in% c("d","p","r","q","s"))
    args$type = type
    fun = eval(parse(text = paste0(args$family, ".likelihood")))
    return(fun(args))
}

inla.likelihood.parser = function(arg_string){
    args_raw = eval(parse(text = arg_string))
    args_res = list()
    if(args_raw$family == "inla.surv"){
        args_raw = args_raw$family.arg.str
    }
    args_res$family = args_raw$family
    args_res$linear.predictor = args_raw$linear.predictor
    args_res$theta = args_raw$theta
    args_res$inv.link.function = eval(parse(text = paste0("inla.link.inv",tolower(args_raw$link.model))))
    args_res$log = FALSE
    args_res$lower.tail = TRUE
    #For Weibulll, we only allow variant 1.
    args_res$variant = args_raw$variant
    #For poisson ... 
    args_res$E = args_raw$E
    #For Gaussian
    args_res$scale = args_raw$scale
    return(args_res)
}







