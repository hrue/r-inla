rMatern <- function(n, coords, kappa, variance, nu=1) {
    m <- as.matrix(dist(coords))
    m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
    diag(m) <- 1
    return(drop(crossprod(chol(variance*m),
                          matrix(rnorm(nrow(coords)*n), ncol=n))))
}

rspde <- function(coords, kappa, variance=1, alpha=2, n=1, mesh, 
                  verbose=FALSE, seed, return.attributes=FALSE) {
    t0 <- Sys.time()
    theta <- c(-0.5*log(4*pi*variance*kappa^2), log(kappa))
    if (verbose) cat('theta =', theta, '\n')
    if (missing(mesh)) {
        mesh.pars <- c(0.5, 1, 0.1, 0.5, 1)*sqrt(alpha-ncol(coords)/2)/kappa 
        if (verbose) cat('mesh.pars =', mesh.pars, '\n')
        attributes <- list(
            mesh=inla.mesh.2d(,
                coords[chull(coords), ], max.edge=mesh.pars[1:2], 
                cutoff=mesh.pars[3], offset=mesh.pars[4:5]))
        if (verbose) cat('n.mesh =', attributes$mesh$n, '\n')
    }
    else attributes <- list(mesh=mesh)
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha=alpha)
    attributes$Q <- inla.spde2.precision(attributes$spde, theta=theta)
    attributes$A <- inla.mesh.project(mesh=attributes$mesh, loc=coords)$A
    if (n==1) 
        result <- drop(attributes$A%*%inla.qsample(
            Q=attributes$Q,
            constr=attributes$spde$f$extraconstr))
    t1 <- Sys.time() 
    result <- inla.qsample(n, attributes$Q, 
                           seed=ifelse(missing(seed), 0, seed), 
                           constr=attributes$spde$f$extraconstr) 
    if (nrow(result)<nrow(attributes$A)) {
        result <- rbind(result, matrix(
            NA, nrow(attributes$A)-nrow(result), ncol(result)))
        dimnames(result)[[1]] <- paste('x', 1:nrow(result), sep='')
        for (j in 1:ncol(result)) 
            result[, j] <- drop(attributes$A%*%
                                result[1:ncol(attributes$A),j])
    }
    else {
        for (j in 1:ncol(result)) 
            result[1:nrow(attributes$A), j] <-
                drop(attributes$A%*%result[,j]) 
        result <- result[1:nrow(attributes$A), ]
    }
    t2 <- Sys.time()
    attributes$cpu <- c(prep=t1-t0, sample=t2-t1, total=t2-t0)
    if (return.attributes) 
        attributes(result) <- c(attributes(result), attributes)
    return(drop(result))
}

precDetFun <- function(U, Q, r) {
  d1 <- determinant(crossprod(U/r, U) + Q)$modulus
  return(new('numeric', (d1 - determinant(Q)$modulus)
             + nrow(U)*log(r))) 
}

precFun <- function(U, Q, r) { 
  Bi <- Diagonal(nrow(U), rep(1/r, nrow(U))) 
  VBi <- crossprod(U, Bi) 
  R <- solve(Q + VBi%*%U, VBi) 
  forceSymmetric(Bi - crossprod(VBi, R)) 
}

negLogLikFun <- function(pars, X, A, y, spde) {
  m <- inla.spde2.precision(spde, c(-pars[2]-1.265512, pars[2])) 
  ldet <- precDetFun(A, m, exp(pars[1])) 
  m <- precFun(A, m, exp(pars[1]))
  Xm <- crossprod(X,m)
  betah <- drop(solve(Xm%*%X, Xm%*%y)) 
  z <- drop(y-X%*%betah)
  s2x.h <- mean(crossprod(m, z)*z) 
  return((ldet + nrow(m)*(1+log(2*pi*s2x.h)))/2) 
}

par2user <- function(pars, X, A, y, spde) {
  m <- inla.spde2.precision(spde, c(-pars[2]-1.265512, pars[2]))
  m <- precFun(A, m, exp(pars[1])) 
  Xm <- crossprod(X, m) 
  beta <- drop(solve(Xm%*%X, Xm%*%y)) 
  z <- drop(y-X%*%beta) 
  s2x.h <- mean(crossprod(m, z)*z)
  c(beta=beta, s2e=exp(pars[1])*s2x.h, 
    s2x=s2x.h, kappa=exp(pars[2])) 
}
