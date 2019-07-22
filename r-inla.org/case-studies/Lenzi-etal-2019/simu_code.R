####-----------------------------------------------------------------####
#### Code for the simulated example in Section 3 of Lenzi et al 2019 ####
####-----------------------------------------------------------------####

# dependencies
library(INLA)
source("auxiliary_functions.R")

# set the input values
sig.eps = 0.5  ### fixed precision for the noise
n.sim = 25    ### number of simulated time series 
N=40          ### length of each time series

s1 = 1:n.sim
phis_100 = sin(s1/(5*pi))^2 ### values of phi's used for simulation

# remove values that are too close to zero or one for stability
phis_100 = ifelse(phis_100 < 0.01, 0.01, phis_100)
phis_100 = ifelse(phis_100 > 0.99, 0.99, phis_100)

# simulating the data
set.seed(1)
y = matrix(c(NA), N, n.sim)

for(ind.sim in 1:n.sim){
  y[,ind.sim] = as.vector(data_y_sim(N=N, phi=phis_100[ind.sim], sig.x=1, sig.eps=sig.eps))
}

##################
##### STEP 1 #####
##################
step1_all = step1(data=y, N=N, sig.x=1, sig.eps=sig.eps, n.sim=n.sim)
# phi's in original scale for plotting
step1_mode_phi = step1_all$mode_phi_first
# EMLCPO
step1_emlcpo = exp(mean(log(unlist(step1_all$cpo_all))))

##################
##### STEP 2 #####
##################
# compare different levels of smoothing on a grid
# par is log(prec) - the higher the smoother
par_all = seq(-5,15,length=5)
step2_all = list(NULL)

for(ind in 1:length(par_all)){
  step2_all[[ind]] = step2(par=par_all[ind], n.sim=n.sim, integration=TRUE,
                           mode_phi_first_intern=step1_all$mode_phi_first_intern, sd_phi_first_intern=step1_all$sd_phi_first_intern)
  print(ind)
}

# plot estimated phi's before and after smoothing
par(mar=c(5,5,4,2)+.1)
plot(s1, phis_100, type="b", xlab="r", ylab=expression(paste(phi)[r]), pch=15, cex.lab=2.5, cex.axis=2.5, lwd=2, ylim=c(-1,1))
points(s1, step1_mode_phi , type="b", lwd=1, col="red", pch=16)
points(s1, step2_all[[2]]$mean_phi_smooth, type="b", lwd=1, col="blue", pch=17)


##################
##### STEP 3 #####
##################
step3_all = list(NULL)
for(ind in 1:length(par_all)){
  step3_all[[ind]] = step3(data=y, mode_phi_smooth_intern=step2_all[[ind]]$mean_phi_smooth_intern, 
                           N=N, sig.x=1, sig.eps=sig.eps, n.sim=n.sim, integration=TRUE, design.data=step2_all[[ind]]$design)
  print(ind)
}

# EMLCPO
step3_emlcpo = rep(NA, length(par_all))
for(ind in 1:length(par_all)){
  step3_emlcpo[ind] =  exp(mean(log(unlist(step3_all[[ind]]))))
  print(ind)
}

####--------------------------------------------------###
#### K-L distance for degrees of smoothness on a grid ### 
#### the lower the better                             ###
####--------------------------------------------------###

kld_step1_all = rep(NA, n.sim)
kld_step3_all = matrix(NA, length(par_all), n.sim)

for(ind in 1:n.sim){
  
  # Sigma_0 and mu_0 for the true (simulated)
  Qx = ar1.prec(N = N, sigma.x = 1, phi = phis_100[ind])
  params_true = mean_sd_sim(ind.sim=ind, data=y, Q=Qx, sig.eps=sig.eps)
  
  # Sigma_1 and mu_1 before smoothing (Step 1)
  step1_simu = step1_mean_cov(data=y, ind.sim=ind, sig.x=1, sig.eps=sig.eps, nn=1000)
  
  # KL-divergence before smoothing 
  kld_step1_all[ind] = KL_distance_multi(mean0_vec=params_true$mean, cov0_mat=params_true$Sigma, 
                                         mean1_vec=step1_simu$mean_mat, cov1_mat=step1_simu$cov_mat)[1]
  
  for(ind.par in 1:length(par_all)){
    
    # Sigma_1 and mu_1 for degrees of smoothness on a grid (Step 3)
    step3_simu = step3_mean_cov(mode_phi_smooth_intern=step2_all[[ind.par]]$mean_phi_smooth_intern, 
                                data=y, ind.sim=ind, sig.x=1, sig.eps=sig.eps, nn=1000) 
    # KL-divergence after smoothing 
    kld_step3_all[ind.par, ind] = KL_distance_multi(mean0_vec=params_true$mean, cov0_mat=params_true$Sigma, 
                                                    mean1_vec=step3_simu$mean_mat, cov1_mat=step3_simu$cov_mat)[1]
    
    rm(step3_simu)
  }
  print(ind)
}

emlkl_step1 = exp(mean(log(kld_step1_all)))
emlkl_step3 = apply(kld_step3_all, 1, function(x) {exp(mean(log(x)))})

####-------------------------------------------------------------------------###
#### Figure 4 of the manuscript: EMLCPO and EMLKL against level of smoothing ###
####-------------------------------------------------------------------------###

# EMLCPO
xx = c(par_all[1]-(par_all[2]-par_all[1]), par_all)
cpos_smooths_diff = c(step1_emlcpo, step3_emlcpo)

inla.dev.new()
par(mar=c(5,5,4,2)+.1)
plot(xx, cpos_smooths_diff, type="b", ylab="EMLCPO", xlab=expression(paste("log(", tau[u], ")")),
     lty=1, pch=15, cex.lab=2, cex.axis=2, lwd=2, xaxt="n")
xx2 = c("no smooth", as.character(round(xx[-1],0)))
axis(1, at=xx, labels=xx2, cex.axis=2)


inla.dev.new()
# EMLKL 
kld_smooths_diff = c(emlkl_step1, emlkl_step3)

par(mar = c(5,5,2,5))
plot(xx, kld_smooths_diff, type="b", ylab="EMLKL", xlab=expression(paste("log(", tau[u], ")")),xlim=range(xx) + c(-0.5, 0.5),
     lty=2, pch=15, cex.lab=2, cex.axis=2, lwd=3, xaxt="n")
xx2 = c("no smooth", as.character(round(xx[-1],0)))
axis(1, at=xx, labels=xx2, cex.axis=2)
