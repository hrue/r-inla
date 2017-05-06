## Export: inla.hyperpar.sample

##!\name{inla.hyperpar.sample}
##!\alias{inla.hyperpar.sample}
##!\alias{inla.hyperpar.sampler}
##!\alias{hyperpar.sample}
##!\alias{hyperpar.sampler}
##!
##!\title{Produce samples from the approximated joint posterior for the hyperparameters}
##!
##!\description{Produce samples from the approximated joint posterior for the hyperparameters}
##!
##!\usage{
##!inla.hyperpar.sample(n, result, intern=FALSE, improve.marginals = FALSE)
##!}
##!
##!\arguments{
##!  \item{n}{Integer. Number of samples required.}
##!  \item{result}{An \code{inla}-object,  f.ex the output from an \code{inla}-call.}
##!  \item{intern}{Logical. If \code{TRUE} then produce samples in the
##!        internal scale for the hyperparmater, if \code{FALSE} then produce
##!        samples in the user-scale. (For example log-precision (intern)
##!        and precision (user-scale))}
##!  \item{improve.marginals}{Logical. If \code{TRUE}, then improve the samples taking into account
##!        possible better marginal estimates for the hyperparameters in \code{result}.} 
##!}
##!
##!\value{%%
##! A matrix where each sample is a row. The contents of the column
##! is described in the rownames.
##!}
##!%%
##!
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!
##!\examples{
##!n = 100
##!r = inla(y ~ 1 + f(idx), data = data.frame(y=rnorm(n), idx = 1:n))
##!ns = 500
##!x = inla.hyperpar.sample(ns, r)
##!
##!rr = inla.hyperpar(r)
##!xx = inla.hyperpar.sample(ns, rr, improve.marginals=TRUE)
##!}


`inla.hyperpar.sample` = function(n, result, intern=FALSE, improve.marginals = FALSE)
{
    ## generate 'n' samples from the joint for the hyperparameters
    ## computed using the CCD approach. if intern==TRUE, then generate
    ## the variables in their intern scale otherwise in the user
    ## scale.
    
    stopifnot(!is.null(result))
    stopifnot(any(class(result) == "inla"))
    stopifnot(n > 0)

    if (improve.marginals) {
        ns = max(n, 300L)
    } else {
        ns = n
    }

    sigma = result$misc$cov.intern
    if (dim(sigma)[1L] == 0L) {
        return (NULL)
    }
    
    p = dim(sigma)[1L]
    sd.plus = result$misc$stdev.corr.positive
    sd.neg  = result$misc$stdev.corr.negative

    if (is.null(sd.plus) || is.null(sd.neg)) {
        sd.plus = rep(1.0, p)
        sd.neg = rep(1.0, p)
    }

    ## do each column  dim(z) = ns x p,  each row is one sample
    z = matrix(NA, ns, p)
    for(i in 1L:p) {
        prob = c(sd.plus[i], sd.neg[i])
        direction = sample.int(2L, ns, prob = prob, replace=TRUE)
        s = c(sd.plus[i], -sd.neg[i])
        z[, i] = s[direction] * abs(rnorm(ns))
    }
    A = result$misc$cov.intern.eigenvectors %*%
        diag(sqrt(result$misc$cov.intern.eigenvalues), nrow = p, ncol = p)
    theta = apply(z, 1L,
                  function(x, A, m) A %*% x + m,
                  A = A, m = result$misc$theta.mode)
    if (p > 1L) {
        theta = t(theta)
    } else {
        theta = matrix(theta, ns, p)
    }
    
    if (improve.marginals) {
        for(i in 1L:p) {
            cdf.old = ecdf(theta[,i]) 
            theta[,i] <- inla.qmarginal((ns/(ns+1.0))* cdf.old(theta[,i]),
                                        result$internal.marginals.hyperpar[[i]])
        }
    }

    theta = theta[1:n,, drop=FALSE]
    if (!intern) {
        for(i in 1L:p) {
            theta[, i] = result$misc$from.theta[[i]](theta[, i])
        }
    }

    if (intern) {
        colnames(theta) = names(result$misc$to.theta)
    } else {
        colnames(theta) = inla.transform.names(result, names(result$misc$to.theta))
    }
    rownames(theta) = paste("sample-", inla.num(1L:n), sep="")

    return (theta)
}

inla.transform.names = function(res, nms)
{
    ## enter name of a hyperpar, return the corresponding name in either in the internal or user
    ## scale
    ret = c()
    for (nm in nms) {
        def = paste(nm, "-- unknown")
        if (is.null(res$internal.summary.hyperpar) || is.null(res$summary.hyperpar)) {
            ret = c(ret, def)
        } else {
            i = which(nm == rownames(res$internal.summary.hyperpar))
            if (length(i) == 1) {
                ret = c(ret, rownames(res$summary.hyperpar)[i])
            } else {
                i = which(nm == rownames(res$summary.hyperpar))
                if (length(i) == 1) {
                    ret = c(ret, rownames(res$internal.summary.hyperpar)[i])
                } else {
                    ret = c(ret, def)
                }
            }
        }
    }
    return (ret)
}
