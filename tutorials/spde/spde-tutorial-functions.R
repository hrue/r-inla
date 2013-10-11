rspde <- function(coords, kappa, variance=1, alpha=2, n=1, 
                  verbose=TRUE, mesh, return.attributes=TRUE) {
    t0 <- Sys.time()
    theta <- c(NA, log(kappa))
    theta[1] <- -0.5*log(4*pi*variance*kappa^2)
    if (verbose) cat('theta =', theta, '\n')
    if (missing(mesh)) {
        mesh.pars <- 1/(kappa*c(3/2,1,7,3,1)) 
        if (verbose) cat('mesh.pars =', mesh.pars, '\n')
        attributes <- list(
            mesh=inla.mesh.2d(
                coords, max.edge=mesh.pars[1:2], 
                cutoff=mesh.pars[3], offset=mesh.pars[4:5]))
        if (verbose) cat('n.mesh =', attributes$mesh$n, '\n')
    }
    else attributes <- list(mesh=mesh)
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha=alpha)
    attributes$Q <- inla.spde2.precision(attributes$spde, theta=theta)
    attributes$A <- inla.mesh.project(mesh=attributes$mesh, loc=coords)$A
    if (n==1) 
        result <- drop(attributes$A%*%inla.qsample(Q=attributes$Q))
    t1 <- Sys.time() 
    result <- inla.qsample(n, attributes$Q)
    if (nrow(result)<nrow(attributes$A)) {
        result <- rbind(result, matrix(
            NA, nrow(attributes$A)-nrow(result), ncol(result)))
        dimnames(result)[[1]] <- paste('x', 1:nrow(result), sep='')
        for (j in 1:ncol(result)) 
            result[, j] <- drop(attributes$A%*%result[1:ncol(attributes$A),j])
    }
    else {
        for (j in 1:ncol(result)) 
            result[1:nrow(attributes$A), j] <- drop(attributes$A%*%result[,j]) 
        result <- result[1:nrow(attributes$A), ]
    }
    t2 <- Sys.time()
    attributes$cpu <- c(prep=t1-t0, sample=t2-t1, total=t2-t0)
    if (return.attributes) 
        attributes(result) <- c(attributes(result), attributes)
    return(result)
}

prec.det.f <- function(U, Q, tau.y) {
  d1 <- determinant(crossprod(U*tau.y, U) + Q)$modulus
  (d1 - determinant(Q)$modulus) - nrow(U)*log(tau.y)
}

precy.f <- function(U, Q, tau.y) {
  Bi <- Diagonal(nrow(U), rep(tau.y, nrow(U)))
  VBi <- crossprod(U, Bi)
  R <- inla.qinv(Q + VBi%*%U) 
  forceSymmetric(Bi - crossprod(VBi, R)%*%VBi)
}

nllf <- function(pars, X, y, A, verbose=TRUE) {
  if (verbose) cat(pars, '')
  m <- inla.spde2.precision(spde, theta=pars[1:2]) 
  ldet <- prec.det.f(A, m, exp(-pars[3]))
  m <- precy.f(A, m, exp(-pars[3]))
  Xm <- crossprod(X,m)
  beta <- drop(solve(Xm%*%X, sum(Xm*y)))
  y <- drop(y-X%*%beta)
  ss <- y%*%m%*%y
  r <- drop((ss + length(y)*log(2*pi) + ldet)/2)
  if (verbose) cat(beta, r, '\n')
  attr(r, 'beta') <- beta
  return(r)
}

par2user <- function(pars, A, X, y) {
  m <- inla.spde2.precision(spde, theta=pars[1:2]) 
  m <- precy.f(A, m, exp(-pars[3]))
  Xm <- crossprod(X,m)
  beta <- drop(solve(Xm%*%X, sum(Xm*y)))
  k2s22 <- exp(sum(pars[1:2])*2)
  c(s2x=1/(4*pi*k2s22), kappa=exp(pars[2]), s2e=exp(pars[3]), beta=beta)
}
