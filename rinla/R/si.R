## Nothing to export. Experimental code....


`inla.si` = function(result, quantiles = c(0.05), n.samples = 10000, debug = FALSE)
{
    if (debug) {
        print("enter inla.si")
    }
    stopifnot(inla.require("mvtnorm"))
    stopifnot(!is.null(result$si))

    quantiles = quantiles[ quantiles < 0.5 ]
    if (length(quantiles) == 0)
        return (result)

    result.new = result
    rc = result$si
    nc = length(rc)
    
    for (itag in 1:length(rc[[1]]$tag)) {
        n = rc[[1]]$len[itag]
        tag= rc[[1]]$tag[itag]
        xx = matrix(0, n, n)
        x = numeric(n)
        x[] = 0

        if (debug) {
            cat("\n")
            cat("\titag = ", itag, "\n")
            cat("\ttag = ", tag, "\n")
            cat("\tn = ", n, "\n")
            cat("\tstart = ", rc[[1]]$start[itag], "\n")
            cat("\tlen = ", rc[[1]]$len[itag], "\n")
        }
        p = numeric(nc)
        for(config in 1:nc)
            p[config] = rc[[config]]$log.dens
        p = exp(p - max(p))
        p = p/sum(p)
        
        for(config in 1:nc) {
            if (debug) {
                cat("\tconfig=", config, "prob=", p[config],"\n")
            }

            m = rc[[config]]$mean[[itag]]
            s = rc[[config]]$sd[[itag]]
            cor = matrix(rc[[config]]$cor[[itag]], n, n)

            xx = xx + p[config] * ( diag(s) %*% cor %*% diag(s) + m %*% t(m) )
            x = x + p[config] * m
        }
        mean = x
        cov  = xx - x %*% t(x)
        s = sqrt(diag(cov))
        cor = diag(1/s) %*% cov %*% diag(1/s)

        for (alpha in sort(quantiles)) {
            sol = uniroot(
                    function(gamma, alpha) {
                        n = dim(cor)[1]
                        lower = rep(qnorm(gamma/2), n)
                        upper = -lower
                        return (pmvnorm(lower = lower, upper = upper, corr = cor)/(1-alpha) - 1)
                    },
                    interval = c(alpha/n, alpha), tol = 1e-3, alpha=alpha)
            gamma = sol$root
            if (debug) {
                cat("\tsolve for alpha=", alpha, "solution gamma = ", gamma, "iter=", sol$iter, "f.root=", sol$f.root,"\n")
            }

            id = which(names(result$marginals.random)== tag)
            stopifnot(length(id) == 1)
            margs = result$marginals.random[[id]]
            nn = length(margs)
            if (debug) {
                cat("\tid=", id, "nn=", nn, "\n")
            }
            lower = numeric(nn)
            upper = numeric(nn)
            for(i in 1:nn) {
                lower[i] = inla.qmarginal(gamma/2, margs[[i]])
                upper[i] = inla.qmarginal(1-gamma/2, margs[[i]])
            }

            id = which(names(result$summary.random)== tag)
            stopifnot(length(id) == 1)
            
            result.new$summary.random[[id]] = cbind(result.new$summary.random[[id]], lower=lower, upper=upper)
            nm = colnames(result.new$summary.random[[id]])
            nm[c(length(nm)-1, length(nm))] = c(paste(1-alpha, "si.low.quant", sep=""),
                      paste(1-alpha, "si.upp.quant", sep=""))
            colnames(result.new$summary.random[[id]]) = nm
            if (debug) {
                cat("\tadded si.interval for alpha=", alpha, "\n")
            }
        }

        ## test the zero-vector
        mm = list()
        cc = list()
        for(config in 1:nc) {
            mm[[config]] = rc[[config]]$mean[[itag]]
            s = rc[[config]]$sd[[itag]]
            cor = matrix(rc[[config]]$cor[[itag]], n, n)
            cc[[config]] = diag(s) %*% cor %*% diag(s)
        }
        ## pvalue = inla.si.mixture(p, mm, cc, rep(0, n), n.samples = n.samples)
    }
    return (result.new)
}
`inla.si.mixture` = function(p, mm, cc, x.star, n.samples, diagonal = 1.0e-4)
{
    n = length(x.star)
    n.mix = length(p)
    nn = as.integer(n.samples * p)
    nn = pmax(1, nn)
    n.samples = sum(nn)

    d = numeric(n.samples)
    j=0
    for(i in 1:n.mix) {
        print(i)
        print(mm[[i]])
        print(cc[[i]])
        print(eigen(cc[[i]], symmetric=TRUE, only.values = TRUE)$values)

        xx = rmvnorm(nn[i], mean = mm[[i]], sigma = cc[[i]] + diag(diagonal, nrow=n))
        dd = numeric(nn[i])
        dd[] = 0
        for(k in 1:n.mix) {
            dd = dd + p[k] * dmvnorm(xx, mean = mm[[k]], sigma = cc[[k]] + diag(diagonal, nrow=n), log=FALSE)
        }
        d[j + 1:nn[i]] = dd
        j = j + nn[i]
    }

    dd = 0
    for(k in 1:n.mix) {
        dd = dd + p[k] * dmvnorm(x.star, mean = mm[[k]], sigma = cc[[k]], log=FALSE)
    }
    stopifnot(length(d) == n.samples)
    p.value = sum(dd < d)/n.samples

    return (p.value)
}
