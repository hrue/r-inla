###--------------------------------------------------###
### AUXILIRARY FUNCTIONS for the AR(1) example       ###    
###--------------------------------------------------###

## Function to create AR1 precision matrix
ar1.prec = function(N , sigma.x, phi) {
  Ui = sparseMatrix(1:N, 1:N, 
                    x = rep(1, N),
                    dims = c(N,N))
  Uim1 = sparseMatrix(2:N, 1:(N-1), 
                      x = rep(1, N-1),
                      dims = c(N,N))
  T1 = Ui - phi*Uim1
  T1[1,1] = sqrt(1-phi^2)
  Q1 = t(T1) %*% T1
  Q1 = Q1*1/(1-phi^2) * 1/sigma.x^2
  return(Q1)
}

## Function to simulate the data
data_y_sim = function(N, phi, sig.x=1, sig.eps){
    
      ## Creating the ar1 precision matrix
      Qx = ar1.prec(N = N, sigma.x = 1, phi = phi)
      ## Creating the iid noise precision matrix
      Qeps = sparseMatrix(1:N, 1:N, x = rep(1, N), dims = c(N,N))
  
      ## x (unknown) is the underlying signal that we want to reconstruct
      Chol.x = chol(Qx)
      x.sim = solve(Chol.x, rnorm(n = N))
      
      ## eps (unknown) is the simulated noise 
      Chol.eps = chol(Qeps)
      eps.sim = solve(Chol.eps, rnorm(n = N))
      
      ## y is the observed data
      y = sig.x*x.sim + sig.eps*eps.sim

  return(y)
}

#############################################
# STEP 1 :: the model fitted to each region #
############################################

step1 = function(data, N, sig.x=1, sig.eps, n.sim){
  
  mode_phi_first_intern <- sd_phi_first_intern <- rep(NA, n.sim)
  mode_phi_first <- rep(NA, n.sim)
  cpo_all = list(NULL)
  
  for(ind.sim in 1:n.sim){
  
    r = inla(y ~ -1 + f(idx, model="ar1", hyper = list(prec = list(initial = log(1/sig.x^2), fixed=TRUE))),
           data = data.frame(y=data[,ind.sim], idx=1:N),
           family = "gaussian",
           control.predictor = list(compute = TRUE),
           num.threads=1,
           control.compute = list(cpo=TRUE, config=TRUE),
           control.family = list(hyper = list(prec = list(initial = log(1/sig.eps^2), fixed=TRUE))))
  
    # mode in the original scale (for plotting)
    mode_phi_first[ind.sim] = r$summary.hyperpar$mode
  
    # mode and sd in the internal scale (for smoothing)
    mode_phi_first_intern[ind.sim] = r$internal.summary.hyperpar$mode
    sd_phi_first_intern[ind.sim] = r$internal.summary.hyperpar$sd
  
    #cpo
    cpo_all[[ind.sim]] = r$cpo$cpo

    print(ind.sim)
  }
  
  return(list(mode_phi_first=mode_phi_first,
              mode_phi_first_intern=mode_phi_first_intern, sd_phi_first_intern=sd_phi_first_intern,
              cpo_all=cpo_all))
}
  
  
#################################
# STEP 2 :: smoothing the phi's #
#################################

step2 = function(par, mode_phi_first_intern, sd_phi_first_intern, n.sim, integration){

  scale = 1/sd_phi_first_intern^2

  r2 = inla(y ~ 1 + f(idx, model="rw2", constr=T),
          data = data.frame(y=mode_phi_first_intern, idx=1:n.sim),
          family = "gaussian",
          control.predictor = list(compute = TRUE),
          scale=scale, 
          control.mode = list(theta=par, fixed=TRUE),
          control.family = list(hyper =  list(prec = list(initial = 0, fixed=TRUE))))

  # smoothed fitted values in the internal scale
  mean_phi_smooth_intern = r2$summary.fitted.values$mean
  sd_phi_smooth_intern = r2$summary.fitted.values$sd

  # smoothed fitted values in the original scale
  mean_phi_smooth = rep(NA, n.sim)
  for(i in 1:n.sim){
    mean_phi_smooth[i] = inla.emarginal(function(x) 2*exp(x)/(1+exp(x))-1, r2$marginals.fitted.values[[i]])
  }
  
  # user-defined integration points: design matrix with nodes and weights 
  nodes.std = c(-1.5, -0.75, 0, 0.75, 1.5) ### nodes for the gauss hermite quadrature (L=5)
  w = c(0.32, 0.75, 1, 0.75, 0.32)        ### weights for the gauss hermite quadrature (L=5)
  design = array(numeric(),c(n.sim,length(nodes.std),2))
  if(integration==TRUE){
    for(j in 1:n.sim){
      design[j,,] = cbind(mean_phi_smooth_intern[j] + sd_phi_smooth_intern[j]*nodes.std, w)
    }
  }

  return(list(mean_phi_smooth_intern=mean_phi_smooth_intern,
              mean_phi_smooth=mean_phi_smooth, design=design))
}
 
  
######################################################################
# STEP 3 :: re-fit the model to each region using the estimated mode #
######################################################################

step3 = function(data, mode_phi_smooth_intern, N, sig.x=1, sig.eps, n.sim, integration, design.data=NULL){
  
  cpo_all = list(NULL)
  
  if(integration==FALSE){
    for(ind.sim in 1:n.sim){
      
      r3 = inla(y ~ -1 + f(idx, model="ar1", hyper = list(prec = list(initial = log(1/sig.x^2), fixed=TRUE))),
                data = data.frame(y=data[,ind.sim], idx=1:N),
                family = "gaussian",
                control.mode = list(theta=mode_phi_smooth_intern[ind.sim], fixed=TRUE),
                control.predictor = list(compute = TRUE),
                control.compute=list(cpo=TRUE, config=TRUE),
                control.family = list(hyper = list(prec = list(initial = log(1/sig.eps^2), fixed=TRUE))))
      
      cpo_all[[ind.sim]] = r3$cpo$cpo
    }
  }
  
  if(integration==TRUE){
    for(ind.sim in 1:n.sim){
      r3 = inla(y ~ -1 + f(idx, model="ar1", hyper = list(prec = list(initial = log(1/sig.x^2), fixed=TRUE))),
                data = data.frame(y=data[,ind.sim], idx=1:N),
                family = "gaussian",
                control.inla = list(int.strategy = "user.expert", int.design = design.data[ind.sim,,]),
                control.predictor = list(compute = TRUE),
                control.compute=list(cpo=TRUE, config=TRUE),
                control.family = list(hyper = list(prec = list(initial = log(1/sig.eps^2), fixed=TRUE))))
      
      cpo_all[[ind.sim]] = r3$cpo$cpo
    }
  }
 
  return(cpo_all)
}



## Compute mean and cov for the true model 

mean_sd_sim = function(ind.sim, data, Q, sig.eps){
  
  tau=1/sig.eps^2
  n=nrow(data)
  
  P = Q + tau*diag(n) 
  b = tau*data[,ind]
  
  mu = solve(P)%*%b
  sigma = solve(P)
  
  return(list(mean = mu, Sigma=sigma))
}



## compute mean and cov for KL-divergence (STEP 1) 

step1_mean_cov = function(data, ind.sim, sig.x, sig.eps, nn){
    
    A = diag(nrow(data))
    lcs = inla.make.lincombs(idx = A)

    r = inla(y ~ -1 + f(idx, model="ar1", hyper = list(prec = list(initial = log(1/sig.x^2), fixed=TRUE))),
             data = data.frame(y=data[,ind.sim], idx=1:nrow(data)),
             family = "gaussian",
             control.predictor = list(compute = TRUE),
             num.threads=1,
             lincomb = lcs, 
             control.inla = list(lincomb.derived.correlation.matrix=TRUE),
             control.compute = list(cpo=TRUE, config=TRUE),
             control.family = list(hyper = list(prec = list(initial = log(1/sig.eps^2), fixed=TRUE))))
    
    m = nrow(data)
    cov_mat = as.matrix(r$misc$lincomb.derived.covariance.matrix, m, m)
    mean_mat = r$summary.lincomb.derived$mean
    
    return(list(cov_mat=cov_mat, mean_mat=mean_mat))
}



## compute mean and cov for KL-divergence (STEP 3) 

step3_mean_cov = function(mode_phi_smooth_intern, data, ind.sim, sig.x=1, sig.eps, nn){

    A = diag(nrow(data))
    lcs = inla.make.lincombs(idx = A)

    r3 = inla(y ~ -1 + f(idx, model="ar1", hyper = list(prec = list(initial = log(1/sig.x^2), fixed=TRUE))),
              data = data.frame(y=data[,ind.sim], idx=1:nrow(data)),
              family = "gaussian",
              control.mode = list(theta=mode_phi_smooth_intern[ind.sim], fixed=TRUE),
              control.predictor = list(compute = TRUE),
              num.threads=1,
              lincomb = lcs, 
              control.inla = list(lincomb.derived.correlation.matrix=TRUE),
              control.compute=list(cpo=TRUE, config=TRUE),
              control.family = list(hyper = list(prec = list(initial = log(1/sig.eps^2), fixed=TRUE))))

    m = nrow(data)
    cov_mat = matrix(r3$misc$lincomb.derived.covariance.matrix, m, m)
    mean_mat = r3$summary.lincomb.derived$mean
    return(list(cov_mat=cov_mat, mean_mat=mean_mat))
}


## compute multivariate K-L Distance from 2 normals 

KL_distance_multi = function(mean0_vec, mean1_vec, cov0_mat, cov1_mat){
    
    L1 = chol(cov0_mat)
    L2 = chol(cov1_mat)
    
    log.num1 = 2*sum(log(diag(L1))) 
    log.num2 = 2*sum(log(diag(L2))) 
    
    d = (1/2) * ( log.num2 - log.num1 - length(mean0_vec) + sum(diag(solve(cov1_mat, cov0_mat))) + 
                  t(mean1_vec - mean0_vec) %*% solve(cov1_mat, mean1_vec - mean0_vec))
    return(d)
}

      
            
      
      
