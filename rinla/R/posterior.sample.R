## Export: inla.posterior.sample

##! \name{inla.sample}
##! \alias{inla.posterior.sample}
##! \alias{posterior.sample}
##! 
##! \title{Generate samples from an approximated posterior of a fitted model}
##! 
##! \description{This function generate samples from an approximated posterior
##!              of a fitted model (an inla-object}
##! \usage{
##!     inla.posterior.sample(n = 1L, result, intern = FALSE, use.improved.mean = TRUE,
##!                           add.names = TRUE, seed = 0L)
##! }
##! 
##! \arguments{
##!   \item{n}{Number of samples.}
##!   \item{result}{The inla-object, ie the output from an \code{inla}-call.
##!       The \code{inla}-object must be created with
##!       \code{control.compute=list(config=TRUE)}.}
##!   \item{use.improved.mean}{Logical. If \code{TRUE} then use the
##!       marginal mean values when constructing samples. If \code{FALSE}
##!       then use the mean in the Gaussian approximations.}
##!  \item{intern}{Logical. If \code{TRUE} then produce samples in the
##!       internal scale for the hyperparmater, if \code{FALSE} then produce
##!       samples in the user-scale. (For example log-precision (intern)
##!       and precision (user-scale))}
##!   \item{add.names}{Logical. If \code{TRUE} then add name for each elements of each
##!       sample. If \code{FALSE}, only add name for the first sample. 
##!       (This save space.)}
##!   \item{seed}{Control the RNG of \code{inla.qsample},
##!       see \code{?inla.qsample} for further information.
##!       If \code{seed=0L} then GMRFLib will set the seed intelligently/at 'random'.
##!       If \code{seed < 0L}  then the saved state of the RNG will be reused if possible, otherwise,
##!       GMRFLib will set the seed intelligently/at 'random'.
##!       If \code{seed > 0L} then this value is used as the seed for the RNG.
##!       If you want reproducible results,  you ALSO need to control the seed for the RNG in R by
##!       controlling the variable \code{.Random.seed} or using the function \code{set.seed},
##!       the example for how this can be done. }
##!}
##!\details{The hyperparameters are sampled from the configurations used to do the
##!       numerical integration,  hence if you want a higher resolution,  you need to
##!       to change the \code{int.stratey} variable and friends. The latent field is
##!       sampled from the Gaussian approximation conditioned on the hyperparameters,
##!       but with a correction for the mean (default).
##!}
##!\value{A list of the samples, where each sample is a list with
##!     names \code{hyperpar} and \code{latent}, and with their marginal
##!     densities in \code{logdens$hyperpar} and \code{logdens$latent}
##!     and the joint density is in \code{logdens$joint}.
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##! 
##!\examples{
##!  r = inla(y ~ 1 ,data = data.frame(y=rnorm(1)), control.compute = list(config=TRUE))
##!  samples = inla.posterior.sample(2,r)
##!
##!  ## reproducible results:
##!  set.seed(1234)
##!  inla.seed = as.integer(runif(1)*.Machine$integer.max)
##!  x = inla.posterior.sample(100, r, seed = inla.seed)
##!  set.seed(1234)
##!  xx = inla.posterior.sample(100, r, seed = inla.seed)
##!  all.equal(x, xx)
##!}


`inla.posterior.sample` = function(n = 1, result, intern = FALSE,
    use.improved.mean = TRUE, add.names = TRUE, seed = 0L)
{
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

    all.samples = rep(list(c()), n)
    i.sample = 1L
    for(k in 1:cs$nconfig) {
        if (n.idx[k] > 0) {
            ## then the latent field
            xx = inla.qsample(n=n.idx[k], Q=cs$config[[k]]$Q,
                    mu = inla.ifelse(use.improved.mean, cs$config[[k]]$improved.mean, cs$config[[k]]$mean), 
                    constr = cs$constr, logdens = TRUE, seed = seed)
            nm = c()
            ld.theta = cs$max.log.posterior + cs$config[[k]]$log.posterior
            for(j in 1:length(cs$contents$tag)) {
                ii = seq(cs$contents$start[j], length = cs$contents$length[j])

                ## if 'tag' is a f() term,  use the ID-names from there
                tag = cs$contents$tag[j]
                random.idx = which(names(result$summary.random) == tag)
                if (length(random.idx) != 1L) {
                    if (cs$contents$length[j] == 1L) {
                        ## this corresponds to fixed effects,  no need to do "(Intercept):1"
                        ## instead of just "(Intercept)"
                        nm = c(nm, cs$contents$tag[j])
                    } else {
                        nm = c(nm,
                               paste(cs$contents$tag[j],
                                     ":",
                                     inla.ifelse(cs$contents$length[j] == 1L, 1, inla.num(1:cs$contents$length[j])), 
                                     sep=""))
                    }
                } else {
                    nm = c(nm,
                           paste(cs$contents$tag[j],
                                 ":",
                                 result$summary.random[[random.idx]]$ID, 
                                 sep=""))
                }
            }
            
            theta = cs$config[[k]]$theta
            log.J = 0.0
            if (!is.null(theta) && !intern) {
                for(j in 1:length(theta)) {
                    theta[j] = do.call(result$misc$from.theta[[j]], args = list(theta[j]))
                }
                names(theta) = inla.transform.names(result, names(theta))
                
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
                rownames(a.sample$latent) = if (add.names || i.sample == 1L) nm else NULL
                all.samples[[i.sample]] = a.sample
                i.sample = i.sample + 1L
            }    
        }
    }
    
    return (all.samples)
}
