library(mvtnorm)
library(tmvtnorm)
source("explore-prior.R")

dat = read.table("SixCity.data", header=TRUE)

ndata = dim(dat)[1]
n = 4
stopifnot(dim(dat)[2] == n+1)

XX = Y = C = list()
for(i in 1:ndata) {
    age = 7:10
    smoking = dat[i, 1]
    x1 = 1
    x2 = age - 9
    x3 = smoking
    x4 = x2 * x3
    y = dat[i, 2:5]

    XX[[i]] = cbind(x1, x2, x3, x4)
    Y[[i]] = y

    lower = c()
    upper = c()
    for(j in 1:n) {
        if (Y[[i]][j] == 0) {
            lower = c(lower, -Inf)
            upper = c(upper, 0)
        } else {
            lower = c(lower, 0)
            upper = c(upper, Inf)
        }
    }
    C[[i]] = list(lower = lower, upper = upper)
}
Data = list(n=n, ndata=ndata, XX=XX, Y=Y, C=C)

log.acc = function(beta.new, beta.old, R.new, R.old, beta.var, 
        perm = FALSE, lambda = 1, Data, parallel = TRUE,
        mc.cores = 2)
{
    if (perm) {
        theta.new = inla.pc.cormat.R2theta(inla.pc.cormat.permute(R.new))
    } else {
        theta.new = inla.pc.cormat.R2theta(R.new)
    }
    theta.old = inla.pc.cormat.R2theta(R.old)
    log.post = inla.pc.cormat.dtheta(theta.new, lambda=lambda, log=TRUE) -
        inla.pc.cormat.dtheta(theta.old, lambda=lambda, log=TRUE) +
            (-1/2 * sum(beta.new^2)/beta.var) -
                (-1/2 * sum(beta.old^2)/beta.var) 
                
    xx = INLA:::inla.mclapply(1:Data$ndata,
            function(k) {
                return(
                    log(c(ptmvnorm(Data$C[[k]]$lower, Data$C[[k]]$upper,
                                   mean = c(Data$XX[[k]] %*% beta.new),
                                   sigma = R.new))) -
                    log(c(ptmvnorm(Data$C[[k]]$lower, Data$C[[k]]$upper,
                                   mean = c(Data$XX[[k]] %*% beta.old),
                                   sigma = R.old))))
            }, mc.cores = mc.cores, parallel = parallel)
    log.post = log.post + sum(unlist(xx))
    if (!perm) {
        return (log.post)
    } else {
        return (list(lacc = log.post,  R.new = inla.pc.cormat.theta2R(theta.new)))
    }
}

log.acc.exch = function(beta.new, beta.old, R.new, R.old, beta.var,
        lambda = 1, Data, mc.cores = 2, parallel = TRUE)
{
    log.post = (prior.rho.exch(R.new[1, 2], lambda = lambda, p = length(beta.new), log=TRUE) -
                prior.rho.exch(R.old[1, 2], lambda = lambda, p = length(beta.new), log=TRUE))

    log.post = log.post + ((-1/2 * sum(beta.new^2)/beta.var) - (-1/2 * sum(beta.old^2)/beta.var))

    xx = INLA:::inla.mclapply(1:Data$ndata,
            function(k) {
                return(
                    log(c(ptmvnorm(Data$C[[k]]$lower, Data$C[[k]]$upper,
                                   mean = c(Data$XX[[k]] %*% beta.new),
                                   sigma = R.new))) -
                    log(c(ptmvnorm(Data$C[[k]]$lower, Data$C[[k]]$upper,
                                   mean = c(Data$XX[[k]] %*% beta.old),
                                   sigma = R.old))))
            }, mc.cores = mc.cores, parallel = parallel)
    log.post = log.post + sum(unlist(xx))
    return (log.post)
}

make.exch.R = function(n, rho)
{
    ## make an exchangable correlation matrix
    R = matrix(rho, n, n)
    diag(R) = 1
    return (R)
}
        

run.mcmc = function (arg.list)
{
    ## arg.list = list(lambda=, exch=)
    
    lambda = arg.list$lambda
    exch = arg.list$exch
    
    ## initial stuff
    parallel = TRUE
    mc.cores = 2L
    n.mcmc = 12000
    s.theta = 0.05
    s.beta = 0.05
    s.rho = 0.03
        
    beta.old = rep(0, n)
    beta.matrix = matrix(NA, n.mcmc, n)
    beta.prior.var = 0.1 ## as in the Chib-paper
        
    if (exch) {
        R.old = make.exch.R(n, 0.001)
        r.matrix = matrix(NA, n.mcmc, 1)
    } else {
        U = matrix(runif(n^2), n, n)
        UU = U %*% t(U)
        diag(UU) = 0
        R.old = diag(n) + UU*0.1
        r.matrix = matrix(NA, n.mcmc, length(inla.pc.cormat.R2theta(R.old)))
    }
    log.post.old = NULL
    accept = numeric(2)
    accept[] = 0
        
    for(iter in 1:n.mcmc) {

        beta.new = beta.old + rnorm(n, sd = s.beta)
        if (exch) {
            rho.new = R.old[1, 2] + rnorm(1, sd = s.rho)
            R.new = make.exch.R(n, rho.new)
            lacc = log.acc.exch(beta.new, beta.old, R.new, R.old,
                    beta.var = beta.prior.var, 
                    lambda = lambda, Data = Data,
                    mc.cores = mc.cores, parallel = parallel)
        } else {
            theta.old = inla.pc.cormat.R2theta(R.old)
            theta.new = theta.old + rnorm(length(theta.old), sd = s.theta)
            R.new = inla.pc.cormat.theta2R(theta.new)
            lacc = log.acc(beta.new, beta.old, R.new, R.old,
                    beta.var = beta.prior.var, 
                    perm=FALSE, lambda = lambda, Data = Data,
                    mc.cores = mc.cores, parallel = parallel)
        }
    
        acc = exp(min(0, lacc))
        if (runif(1) < acc) {
            accept[1] = accept[1] + 1
            beta.old = beta.new
            R.old = R.new
        }

        ## permute unless exch
        if (!exch) {
            res = log.acc(beta.old, beta.old, R.old, R.old,
                    beta.var = beta.prior.var, 
                    perm=TRUE, lambda = lambda, Data = Data,
                    mc.cores = mc.cores, parallel = parallel)
            if (runif(1) < exp(min(0, res$lacc))) {
                accept[2] = accept[2] + 1
                R.new = R.old = res$R.new
            }
        }

        beta.matrix[iter, ] = c(beta.new)
        if (exch) {
            r.matrix[iter, ] = R.new[1, 2]
        } else {
            r.matrix[iter, ] = inla.pc.cormat.R2r(R.new)
        }

        if (iter %% 10 == 0) {
            print(round(c(iter, accept[1]/iter, accept[2]/iter), 2))
        }
    }

    result = list(
            lambda = lambda,
            exch = exch, 
            beta.matrix = beta.matrix,
            r.matrix = r.matrix,
            accept = accept)
    save(result,
         file = paste("result--lambda-", round(lambda, 4),
                 "--exch-", as.numeric(exch), ".Rsave", sep=""))
}


INLA:::inla.mclapply(list(list(lambda = 0.1, exch = FALSE),
                          list(lambda = 0.1, exch = TRUE),
                          list(lambda = 1, exch = FALSE),
                          list(lambda = 1, exch = TRUE)),
                     run.mcmc, mc.cores = 4)

        
