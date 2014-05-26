## Export: inla.posterior.sample

##! \name{inla.sample}
##! \alias{inla.posterior.sample}
##! \alias{posterior.sample}
##! 
##! \title{Generate samples from an approximated posterior of a fitted model}
##! 
##! \description{This function generate samples from an approximated posterior of a fitted model,  ie an inla-object}
##! \usage{
##!     inla.posterior.sample(n = 1L, result, hyper.user.scale = TRUE, use.improved.mean = TRUE)
##! }
##! 
##! \arguments{
##!   \item{n}{Number of samples.}
##!   \item{result}{The inla-object, ie the output from an \code{inla}-call. The \code{inla}-object must be created with
##!                 \code{control.compute=list(config=TRUE)}.}
##!   \item{hyper.user.scale}{Logical. If \code{TRUE} then values of
##!   the hyperparameters are given in the user scale (for example
##!   \code{precision}). If \code{FALSE} then values of the
##!   hyperparameters are given in the internal representation (for
##!   example \code{log(precision)}).}
##!   \item{use.improved.mean}{Logical. If \code{TRUE} then use the
##!   marginal mean values when constructing samples. If \code{FALSE}
##!   then use the mean in the Gaussian approximations.}
##!}
##!\value{ A list of the samples, where each sample is a list with
##!  names \code{hyperpar} and \code{latent}, and with their marginal
##!  densities in \code{logdens$hyperpar} and \code{logdens$latent}
##!  and the joint density is in \code{logdens$joint}.  THIS IS AN
##!  EXPERIMENTAL FUNCTION AND CHANGES MAY APPEAR AT ANY TIME!  }
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##!  r = inla(y ~ 1 ,data = data.frame(y=rnorm(1)), control.compute = list(config=TRUE))
##!  samples = inla.posterior.sample(2,r)
##!}


`inla.posterior.sample` = function(n = 1, result, hyper.user.scale = TRUE, use.improved.mean = TRUE)
{
    warning("inla.posterior.sample: THIS FUNCTION IS EXPERIMENTAL!!!")
    
    stopifnot(!missing(result) && any(class(result) == "inla"))
    if (is.null(result$misc$configs)) {
        stop("You need an inla-object computed with option 'control.compute=list(config = TRUE)'.")
    }
    n = as.integer(n)
    stopifnot(is.integer(n) && n > 0L)
    
    ## first sample the theta(-index)
    cs = result$misc$configs
    ld = numeric(cs$nconfig)
    for (i in 1:cs$nconfig) {
        ld[i] = cs$config[[i]]$log.posterior
    }
    p = exp(ld - max(ld))
    idx = sample(1:cs$nconfig, n, prob = p, replace = TRUE)
    idx = sort(idx)
    n.idx = numeric(cs$nconfig)
    n.idx[] = 0
    for(i in 1:cs$nconfig) {
        n.idx[i] = sum(idx == i)
    }

    all.samples = list(list())
    all.samples[[n]] = NA
    i.sample = 1L
    for(k in 1:cs$nconfig) {
        if (n.idx[k] > 0) {
            ## then the latent field
            xx = inla.qsample(n=n.idx[k], Q=cs$config[[k]]$Q,
                    mu = inla.ifelse(use.improved.mean, cs$config[[k]]$improved.mean, cs$config[[k]]$mean), 
                    constr = cs$constr, logdens = TRUE)
            nm = c()
            ld.theta = cs$max.log.posterior + cs$config[[k]]$log.posterior
            for(j in 1:length(cs$contents$tag)) {
                ii = seq(cs$contents$start[j], length = cs$contents$length[j])
                nm = c(nm,
                        paste(cs$contents$tag[j],
                              ".",
                              inla.ifelse(cs$contents$length[j] == 1L, 1, inla.num(1:cs$contents$length[j])), 
                              sep=""))
            }
            
            theta = cs$config[[k]]$theta
            log.J = 0.0
            if (!is.null(theta) && hyper.user.scale) {
                for(j in 1:length(theta)) {
                    theta[j] = do.call(result$misc$from.theta[[j]], args = list(theta[j]))
                }
                names(theta) = paste(names(theta), "-- in user scale")
                
                if (TRUE) {
                    ## new fancy code using the automatic differentiation feature in R
                    for(i in 1:length(theta)) {
                        arg.val = formals(result$misc$from.theta[[i]])
                        arg = names(arg.val)
                        if (length(arg) == 1L) {
                            deriv.func = inla.eval(paste("function(", arg, ") {}"))
                        } else {
                            if (length(arg)==2L) {
                                deriv.func = inla.eval(paste("function(", arg[1L], ",",  arg[2L], "=", arg.val[2L], ") {}"))
                            }
                            else {
                                stopifnot(length(arg) == 3L)
                                deriv.func = inla.eval(paste("function(", arg[1L], ",",  arg[2L], "=", arg.val[2L], ",", arg[3L], "=", arg.val[3L], ") {}"))
                            }
                        }
                        temptest <- try(D(body(result$misc$from.theta[[i]]), arg[1L]), silent=TRUE)
                        if (!inherits(temptest, "try-error")) {
                            body(deriv.func) <- temptest 
                            log.J = log.J - log(abs(deriv.func(cs$config[[k]]$theta[i]))) ## Yes, it's a minus...
                        }
                        else {
                            h = .Machine$double.eps^0.25
                            theta.1 = do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] - h))
                            theta.2 = do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] + h))
                            log.J = log.J - log(abs((theta.2 - theta.1)/(2.0*h))) ## Yes, it's a minus...
                        }
                    }
                    ## print(paste("logJ", log.J))
                } else {
                    ## old code using numerical differentiation
                    h = .Machine$double.eps^0.25
                    for(i in 1:length(theta)) {
                        theta.1 = do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] - h))
                        theta.2 = do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] + h))
                        log.J = log.J - log(abs((theta.2 - theta.1)/(2.0*h))) ## Yes, it's a minus...
                    }
                    ## print(paste("logJ", log.J))
                }
            }
            
            for(i in 1:n.idx[k]) {
                if (is.null(theta)) {
                    a.sample = list(
                            hyperpar = NULL,
                            latent = xx$sample[ , i, drop=FALSE],
                            logdens = list(
                                    hyperpar = NULL,
                                    latent = as.numeric(xx$logdens[i]),
                                    joint = as.numeric(xx$logdens[i])))
                } else {
                    ld.h = as.numeric(ld.theta - result$mlik[1, 1] + log.J)
                    a.sample = list(
                            hyperpar = theta,
                            latent = xx$sample[ , i, drop=FALSE],
                            logdens = list(
                                    hyperpar = ld.h,
                                    latent = as.numeric(xx$logdens[i]),
                                    joint = as.numeric(ld.h + xx$logdens[i])))
                }
                rownames(a.sample$latent) = nm
                all.samples[[i.sample]] = a.sample
                i.sample = i.sample + 1L
            }    
        }
    }
    
    return (all.samples)
}
