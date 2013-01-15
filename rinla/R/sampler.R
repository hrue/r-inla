
inla.hyperpar.sampler = function(...) {
    warning("New name is 'inla.hyperpar.sample'")
    return (inla.hyperpar.sample(...))
}

inla.hyperpar.sample = function(n, result, intern=FALSE)
{
    ## generate 'n' samples from the joint for the hyperparameters
    ## computed using the CCD approach. if intern==TRUE, then generate
    ## the variables in their internal scale otherwise in the user
    ## scale.
    
    stopifnot(!is.null(result))
    stopifnot(any(class(result) == "inla"))
    stopifnot(n > 0)

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

    ## do each column  dim(z) = n x p,  each row is one sample
    z = matrix(NA, n, p)
    for(i in 1L:p) {
        ## this is ok as n >> p, usually.
        prob = c(sd.plus[i], sd.neg[i])
        direction = sample.int(2L, n, prob = prob, replace=TRUE)
        s = c(sd.plus[i], -sd.neg[i])
        z[, i] = s[direction] * abs(rnorm(n))
    }
    A = result$misc$cov.intern.eigenvectors %*%
        diag(sqrt(result$misc$cov.intern.eigenvalues), nrow = p, ncol = p)
    theta = apply(z, 1L,
            function(x, A, m) A %*% x + m,
            A = A, m = result$misc$theta.mode)
    ## fix for p > 1
    if (p > 1L) {
        theta = t(theta)
    } else {
        theta = matrix(theta, n, p)
    }
    
    if (!intern) {
        ## map to user-scale. do each column at the time as n >> p,
        ## usually.
        for(i in 1L:p) {
            theta[, i] = result$misc$from.theta[[i]]( theta[, i] )
        }
    }

    if (intern) {
        colnames(theta) = names(result$misc$to.theta)
    } else {
        colnames(theta) = paste(names(result$misc$to.theta), "in user-scale")
    }
    rownames(theta) = paste("sample-", inla.num(1L:n), sep="")

    return (theta)
}
