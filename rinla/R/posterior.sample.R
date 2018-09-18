## Export: inla.posterior.sample inla.posterior.sample.eval

##! \name{inla.sample}
##! \alias{inla.posterior.sample}
##! \alias{posterior.sample}
##! \alias{inla.posterior.sample.eval}
##! \alias{posterior.sample.eval}
##! 
##! \title{Generate samples, and functions thereof, from an approximated posterior of a fitted model}
##! 
##! \description{This function generate samples, and functions of those,
##!              from an approximated posterior of a fitted model (an inla-object)}
##! \usage{
##!     inla.posterior.sample(n = 1L, result, intern = FALSE, use.improved.mean = TRUE,
##!                           add.names = TRUE, seed = 0L, num.threads = 1L)
##!     inla.posterior.sample.eval(fun, samples, return.matrix = TRUE, ...)
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
##!   \item{num.threads}{The number of threads that can be used. \code{num.threads>1L} requires
##!       \code{seed = 0L}. Only use \code{num.threads > 1L} for large problems/number of
##!       samples. This option does currently NOT use the default one set by \code{inla.setOption()}. }
##!   \item{fun}{The function to evaluate for each sample. Upon entry, the variable names
##!               defined in the model are defined as the value of the sample.
##!               The list of names are defined in \code{result$misc$configs$contents} where
##!               \code{result} is an \code{inla}-object. This includes predefined names for
##!               for the linear predictor (\code{Predictor} and \code{APredictor}),  and the
##!               intercept (\code{(Intercept)} or \code{Intercept}).
##!               The hyperparameters are defined as \code{theta},  no matter if they are in the
##!               internal scale or not. The function \code{fun} can also return a vector.}
##!   \item{samples}{\code{samples} is the output from \code{inla.posterior.sample()}}
##!   \item{return.matrix}{Logical. If \code{TRUE},  then return the samples of \code{fun}
##!                         as matrix,  otherwise,  as a list.}
##!   \item{...}{Additional arguments to \code{fun}}
##!}
##!\details{The hyperparameters are sampled from the configurations used to do the
##!       numerical integration,  hence if you want a higher resolution,  you need to
##!       to change the \code{int.stratey} variable and friends. The latent field is
##!       sampled from the Gaussian approximation conditioned on the hyperparameters,
##!       but with a correction for the mean (default).
##!
##!       Set sparse-matrix library with \code{inla.setOption(smtp=...)} and
##!       number of threads by \code{inla.setOption(num.threads=...)}.
##!}
##!\value{\code{inla.posterior.sample} returns a list of the samples,
##!       where each sample is a list with
##!     names \code{hyperpar} and \code{latent}, and with their marginal
##!     densities in \code{logdens$hyperpar} and \code{logdens$latent}
##!     and the joint density is in \code{logdens$joint}.
##!     \code{inla.posterior.sample.eval} return a list or a matrix of 
##!     \code{fun} applied to each sample.
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
##!
##! set.seed(1234)
##! n = 25
##! xx = rnorm(n)
##! yy = rev(xx)
##! z = runif(n)
##! y = rnorm(n)
##! r = inla(y ~ 1 + z + f(xx) + f(yy, copy="xx"),
##!         data = data.frame(y, z, xx, yy), 
##!         control.compute = list(config=TRUE),
##!         family = "gaussian")
##! r.samples = inla.posterior.sample(100, r)
##! 
##! fun = function(...) {
##!     mean(xx) - mean(yy)
##! }
##! f1 = inla.posterior.sample.eval(fun, r.samples)
##! 
##! fun = function(...) {
##!     c(exp(Intercept), exp(Intercept + z))
##! }
##! f2 = inla.posterior.sample.eval(fun, r.samples)
##! 
##! fun = function(...) {
##!     return (theta[1]/(theta[1] + theta[2]))
##! }
##! f3 = inla.posterior.sample.eval(fun, r.samples)
##!
##! ## Predicting nz new observations, and
##! ## comparing the estimated one with the true one
##! set.seed(1234)
##! n = 100
##! alpha = beta = s = 1
##! z = rnorm(n)
##! y = alpha + beta * z + rnorm(n, sd = s)
##! r = inla(y ~ 1 + z, 
##!         data = data.frame(y, z), 
##!         control.compute = list(config=TRUE),
##!         family = "gaussian")
##! r.samples = inla.posterior.sample(10^3, r)
##! nz = 3
##! znew = rnorm(nz)
##! fun = function(zz = NA) {
##!     ## theta[1] is the precision
##!     return (Intercept + z * zz +
##!             rnorm(length(zz), sd = sqrt(1/theta[1])))
##! }
##! par(mfrow=c(1, nz))
##! f1 = inla.posterior.sample.eval(fun, r.samples, zz = znew)
##! for(i in 1:nz) {
##!     hist(f1[i, ], n = 100, prob = TRUE)
##!     m = alpha + beta * znew[i]
##!     xx = seq(m-4*s, m+4*s, by = s/100)
##!     lines(xx, dnorm(xx, mean=m, sd = s), lwd=2)
##! }
##!}

`inla.posterior.sample` = function(n = 1, result, intern = FALSE,
    use.improved.mean = TRUE, add.names = TRUE, seed = 0L, num.threads = 1L)
{
    stopifnot(!missing(result) && any(class(result) == "inla"))
    if (is.null(result$misc$configs)) {
        stop("You need an inla-object computed with option 'control.compute=list(config = TRUE)'.")
    }

    if (is.null(num.threads)) {
        num.threads = inla.getOption("num.threads")
    }
    num.threads = max(num.threads, 1L)
    if (num.threads > 1L) {
        if (seed != 0L) {
            stop("num.threads > 1L require seed = 0L")
        }
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

    con = cs$contents
    all.samples = rep(list(c()), n)
    i.sample = 1L
    for(k in 1:cs$nconfig) {
        if (n.idx[k] > 0) {
            ## then the latent field
            xx = inla.qsample(n=n.idx[k], Q=cs$config[[k]]$Q,
                              mu = inla.ifelse(use.improved.mean, cs$config[[k]]$improved.mean, cs$config[[k]]$mean), 
                              constr = cs$constr, logdens = TRUE, seed = seed, num.threads = num.threads)
            ## if user set seed,  then just continue this rng-stream
            if (seed > 0L) seed = -1L

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
                if (add.names || i.sample == 1L) {
                    n1 = length(nm)
                    n2 = length(a.sample$latent)
                    stopifnot(n2 >= n1)  ## this must be true. just a check
                    if (n2 > n1 ) {
                        ## This is the case where lincomb.derived.only = FALSE, so these are
                        ## then added to the end. Should transfer the names of them all the way
                        ## to here, but...
                        xnm = paste0("Lincomb:", inla.num(1:(n2-n1)))
                        nm = c(nm,  xnm)

                        if (i.sample == 1L) {
                            ## add to the contents if needed
                            con$tag = c(con$tag, "Lincomb")
                            con$start = c(con$start, n1 + 1)
                            con$length = c(con$length, n2 - n1)
                        }
                    }
                    rownames(a.sample$latent) = nm
                } else {
                    rownames(a.sample$latent) = NULL
                }
                all.samples[[i.sample]] = a.sample
                i.sample = i.sample + 1L
            }    
        }
    }

    attr(all.samples, ".contents") = con
    return (all.samples)
}

`inla.posterior.sample.eval` = function(fun, samples, return.matrix = TRUE, ...)
{
    ## evaluate FUN(...) over each sample
    var = "inla.posterior.sample.eval.warning.given"
    if (!(exists(var, envir = inla.get.inlaEnv()) &&
          get(var, envir = inla.get.inlaEnv()) == TRUE)) {
        warning("Function 'inla.posterior.sample.eval()' is experimental.")
        assign(var, TRUE, envir = inla.get.inlaEnv())
    }

    contents = attr(samples, which = ".contents", exact = TRUE)
    if (is.null(contents)) {
        stop("Argument 'samples' must be the output from 'inla.posterior.sample()'.")
    }

    my.fun = function(a.sample, .contents, .fun, ...) 
    {
        env = new.env()
        theta = as.vector(a.sample$hyperpar)
        assign("theta", theta, envir = env)
        for (i in seq_along(.contents$tag)) {
            assign(.contents$tag[i],
                   as.vector(a.sample$latent[.contents$start[i]:(.contents$start[i] +
                                                                 .contents$length[i] - 1), 1]), 
                   envir = env)
        }
        ## this is special
        assign("Intercept", get("(Intercept)", envir = env), envir = env)

        parent.env(env) = .GlobalEnv
        environment(.fun) = env
        return (.fun(...))
    }

    ret = inla.mclapply(samples, my.fun, .fun=fun, .contents = contents, ...)
    if (return.matrix) {
        ns = length(ret)
        ret = matrix(unlist(ret), ncol = ns)
        colnames(ret) = paste0("sample", 1:ns)
        rownames(ret) = paste0("fun", 1:nrow(ret))
    }

    return (ret)
}
