##! \name{inla.sample}
##! \alias{inla.posterior.sample}
##! \alias{posterior.sample}
##! 
##! \title{Generate samples from an approximated posterior of a fitted model}
##! 
##! \description{This function generate samples from an approximated posterior of a fitted model,  ie an inla-object}
##! \usage{
##!     inla.posterior.sample(n = 1L, result, hyper.user.scale = TRUE)
##! }
##! 
##! \arguments{
##!   \item{n}{Number of samples.}
##!   \item{result}{The inla-object, ie the output from an \code{inla}-call. The \code{inla}-object must be created with
##!                 \code{control.compute=list(config=TRUE)}.}
##!   \item{hyper.user.scale}{Use values of hyperparameters in the user scale, not in the internal representation}
##!}
##!\value{
##!  A list of the samples,  where each sample is a list with names \code{hyperpar} and \code{latent},  and with
##!  their marginal densities in \code{logdens$hyperpar} and \code{logdens$latent}
##!  and the joint density is in \code{logdens$joint}.
##!  THIS IS AN EXPERIMENTAL FUNCTION AND CHANGES MAY APPEAR AT ANY TIME!
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##!  r = inla(y ~ 1 ,data = data.frame(y=rnorm(1)), control.compute = list(config=TRUE))
##!  samples = inla.posterior.sample(2,r)
##!}


`inla.posterior.sample` = function(n = 1, result, hyper.user.scale = TRUE)
{
    warning("THIS FUNCTION IS EXPERIMENTAL!!!")
    
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
            tmp = inla.qsample(n=n.idx[k], Q=cs$config[[k]]$Q, constr = cs$constr, logdens = TRUE)
            xx = cs$config[[k]]$mean + tmp$sample
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
            if (hyper.user.scale) {
                for(j in 1:length(theta)) {
                    theta[j] = do.call(result$misc$from.theta[[j]], args = list(theta[j]))
                }
                names(theta) = paste(names(theta), "-- in user scale")
            } 

            for(i in 1:n.idx[k]) {
                a.sample = list(hyperpar = theta,  latent = xx[ , i, drop=FALSE],
                        logdens = list(hyperpar = ld.theta, latent = as.numeric(tmp$logdens[i]),
                                joint = as.numeric(ld.theta + tmp$logdens[i])))
                rownames(a.sample$latent) = nm
                all.samples[[i.sample]] = a.sample
                i.sample = i.sample + 1L
            }    
        }
    }
    
    return (all.samples)
}
