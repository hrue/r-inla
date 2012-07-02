rm(list = ls())

##### loading the needed libraries ################

library (INLA)
library(lme4)
###---------------------------------------------------

###############################################################
###############     Simulation Example     ####################

### Simulating observations from a mixed Poisson count model

set.seed(1234)
simfun <- function(ng = 50, nr = 2, fsd = 0.75, b = c(1,2,-3)) {
	ntot <- nr * ng
	b.reff <- rnorm(ng, sd = fsd)
	x1 <- rbinom(ntot, 1, 0.6) 
	x2 <- runif(ntot, 0, 1) 
	dd <- data.frame(cbind(x1,x2), f = factor(rep(1:ng, each = nr)))
	dd$eta0 <- model.matrix(~cbind(x1,x2), data = dd) %*% b
	dd$eta <- with(dd, eta0 + b.reff[f] )
	dd$mu <- exp(dd$eta)
	dd$y <- with(dd, rpois(ntot, lambda = mu))
	dd
}

dd <- simfun()

###---------------------------------------------------

##### The main simulation function to compare AGHQ and HDC methods #######

cfun.hdc <- function(d) {
#### 
beta0<-matrix(0,ncol=4,nrow=iter)
beta1<-matrix(0,ncol=4,nrow=iter)
beta2<-matrix(0,ncol=4,nrow=iter)
vc<-matrix(0,ncol=4,nrow=iter)
se.beta0<-matrix(0,ncol=4,nrow=iter)
se.beta1<-matrix(0,ncol=4,nrow=iter)
se.beta2<-matrix(0,ncol=4,nrow=iter)
se.vc<-matrix(0,ncol=4,nrow=iter)
###
for(i in 1:iter){
	d<-simfun()
####
##### MLE fitting by adaptive Gauss-Hermit quadrature method with 15 nodes
####
	m <- glmer(y ~ x1 +x2 + (1 | f) , family = "poisson" , data = d , nAGQ=15)
####
##### preparing cloned data vector
####
	"x1.k"<-rep(d[,1],k)
	"x2.k"<-rep(d[,2],k)
	"f.k"<-rep(d[,3],k)
	"y.k"<-rep(d[,7],k)
	clone.data<-cbind(x1.k,x2.k,f.k,y.k)
	clone.data[,3]<-rep(1:(ng*k),each=nr)
	clone.data<-as.data.frame(clone.data)
####
##### HDC algorithm with Informative (gamma) prior for variance component
####
	hdc.1 <- inla(y.k ~ x1.k + x2.k +
	f(f.k,model="iid",param=c(1, 2)),family = "poisson",data = clone.data)
	hyp.hdc.1 = inla.hyperpar(hdc.1)
	vc.mean.1 = inla.expectation(function(x) 1/x^.5, hyp.hdc.1$marginals[[1]])
	vc.m2.1 = inla.expectation(function(x) 1/x, hyp.hdc.1$marginals[[1]])
	stdev.1 = sqrt(vc.m2.1- vc.mean.1^2)

####
##### HDC algorithm with Non-informative (flat) priors for variance component
####
	hdc.2 <- inla(y.k ~ x1.k + x2.k +
	f(f.k,model="iid"),family = "poisson",data = clone.data,
	control.family=list(prior="flat"))
	hyp.hdc.2 = inla.hyperpar(hdc.2)
	vc.mean.2 = inla.expectation(function(x) 1/x^.5, hyp.hdc.2$marginals[[1]])
	vc.m2.2 = inla.expectation(function(x) 1/x, hyp.hdc.2$marginals[[1]])
	stdev.2 = sqrt(vc.m2.2- vc.mean.2^2)
####
##### HDC algorithm with Vague (gamma) priors for variance component
####
	hdc.3 <- inla(y.k ~ x1.k + x2.k +
	f(f.k,model="iid",param=c(3, 0.05)),family = "poisson", data = clone.data)
	hyp.hdc.3 = inla.hyperpar(hdc.3)
	vc.mean.3 = inla.expectation(function(x) 1/x^.5, hyp.hdc.3$marginals[[1]])
	vc.m2.3 = inla.expectation(function(x) 1/x, hyp.hdc.3$marginals[[1]])
	stdev.3 = sqrt(vc.m2.3- vc.mean.3^2)
####
##### Preparing output (parameters)
####
	beta0[i,]<-c(fixef(m)[1],hdc.1$summary.fixed[1,1],hdc.2$summary.fixed[1,1],
	hdc.3$summary.fixed[1,1])
##
	beta1[i,]<-c(fixef(m)[2],hdc.1$summary.fixed[2,1],hdc.2$summary.fixed[2,1],
	hdc.3$summary.fixed[2,1])
##
	beta2[i,]<-c(fixef(m)[3],hdc.1$summary.fixed[3,1],hdc.2$summary.fixed[3,1],
	hdc.3$summary.fixed[3,1])
##
	vc[i,]<-c(sqrt(unlist(VarCorr(m))),vc.mean.1,vc.mean.2,vc.mean.3)
###
###### (precision of parameters)
###
	se.beta0[i,]<-c(sqrt(diag(vcov(m)))[1],hdc.1$summary.fixed[1,2],
	hdc.2$summary.fixed[1,2],hdc.3$summary.fixed[1,2])
##
	se.beta1[i,]<-c(sqrt(diag(vcov(m)))[2],hdc.1$summary.fixed[2,2],
	hdc.2$summary.fixed[2,2],hdc.3$summary.fixed[2,2])
##
	se.beta2[i,]<-c(sqrt(diag(vcov(m)))[3],hdc.1$summary.fixed[3,2],
	hdc.2$summary.fixed[3,2],hdc.3$summary.fixed[3,2])
##
	se.vc[i,]<-c(sqrt(diag(VarCorr(m)$f)/ng),stdev.1,stdev.2,stdev.3)
####
	cat("iter=", i,"\n")
####
	}
	return(list('Beta0'=beta0,'Beta1'=beta1,'Beta2'=beta2,'Sigma'=vc,
	'SE.Beta0'=se.beta0,'SE.Beta1'=se.beta1,'SE.Beta2'=se.beta2,'SE.Sigma'=se.vc))
}

###---------------------------------------------------
####### Run simulation for 100 data sets and 80 clones of data ############
ng = 50; nr = 2
# k=80
# iter=100
# rr.hdc <- cfun.hdc()

###########
###---------------------------------------------------
####### Comparing computing times for DC and HDC methods on  ######################
####### a typical simulated data set 		    		 ######################

k=100

"x1.k"<-rep(dd[,1],k)
"x2.k"<-rep(dd[,2],k)
"f.k"<-rep(dd[,3],k)
"y.k"<-rep(dd[,7],k)
clone.data<-cbind(x1.k,x2.k,f.k,y.k)
clone.data[,3]<-rep(1:(ng*k),each=nr)
clone.data<-as.data.frame(clone.data)


### HDC 

hdc<-inla(y.k ~ x1.k + x2.k + f(f.k,model="iid"), data=clone.data, family="poisson")

### DC
wd.dc = tempfile()
## same as before
result.dc = system.time(inla(y.k ~ x1.k + x2.k + f(f.k,model="iid"),
	data=clone.data, family="poisson",
	working.directory = wd.dc,
  	keep = TRUE,
 	inla.arg = "-m mcmc -N 10000 -T 10 -S 1"))

## Compare cpu times
time.hdc = hdc$cpu.used
time.hdc
result.dc

# Plotting the results  (for example marginals for the hyperparameters)
hyp.hdc = hdc$marginals.hyperpar[[1]]
hyp.dc = scan(paste(wd.dc,"/results.files/hyperparameter-random.effect00000001-parameter-user-scale/trace.dat",sep=""))
plot(hyp.hdc, lty="l", col=1,
main=expression(paste("hybrid DC and DC based densities of",~ sigma^{-2})))
lines(density(hyp.dc), lty="2", col=2)


##### compare another DC and HDC-based densities (for example the slope of x1)

int.inla = hdc$marginals.fixed[[2]]
int.dc = scan(paste(wd.dc,"/results.files/fixed.effect00000002/trace.dat",sep=""))
plot(int.inla, lty="l", col=1,
main=expression(paste("hybrid DC and DC based densities of",~ beta[1])))
lines(density(int.dc), lty="2", col=2)

####--------------------------------------------------------
###################################################################
###############     Overdispersion Example: Seeds Data    #########

################# Preparing cloned data set with k=200 ############

seed.data<-Seeds

n1=21
k=200
n=n1*k

"r.k"<-rep(seed.data[,1],k)
"n.k"<-rep(seed.data[,2],k)
"x1.k"<-rep(seed.data[,3],k)
"x2.k"<-rep(seed.data[,4],k)
"plate.k"<-rep(seed.data[,5],k)

clone.data<-cbind(r.k,n.k,x1.k,x2.k,plate.k)

clone.data[,5]<-1:n
clone.data<-as.data.frame(clone.data)

####---------------------------------------------------------------------------
###################### main effects model  ############################

### Fitting by AGHQ method with 15 nodes #####

(main.effects <- glmer(r/n ~ x1 + x2 + (1 | plate) , 
			family = "binomial", data = Seeds, nAGQ=15))

####---------------------------------------------------------------------------
#### Fitting by INLA method ####

seeds.inla.fit.1 = inla(r ~ x1 + x2 + f(plate, model="iid",
        param=c(.5, .0164)), data=Seeds, family="binomial", Ntrials=n )

seeds.hyperpar.1 = inla.hyperpar(seeds.inla.fit.1)

#### Fixed effects summaries
summary(seeds.inla.fit.1)

#### variance component summaries
mean.1=m1 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.1$marginals[[1]])
m2 = inla.expectation(function(x) 1/x, seeds.hyperpar.1$marginals[[1]])
stdev.1 = sqrt(m2- mean.1^2)

####---------------------------------------------------------------------------
#### Fitting by HDC method with informative prior ######

formula.clone.1 = r.k ~ x1.k+x2.k+f(plate.k,model="iid",param=c(.5,.0164))
mod.seeds.clone.1 = inla(formula.clone.1,data=clone.data,family="binomial",Ntrials=n.k)

seeds.hyperpar.clone.1 = inla.hyperpar(mod.seeds.clone.1)

#### Fixed effects summaries
summary(mod.seeds.clone.1)
sd.clone.1<-as.vector(sqrt(k)*mod.seeds.clone.1$summary.fixed[,2])

#### Preparing elements for plotting 
marginal.clone.intercept<-mod.seeds.clone.1$marginals.fixed[[1]]
marginal.clone.x1<-mod.seeds.clone.1$marginals.fixed$x1
marginal.clone.x2<-mod.seeds.clone.1$marginals.fixed$x2
hyp.k.1<-seeds.hyperpar.clone.1$marginals$`Precision for plate.k`

#### Variance component summaries
mean.clone.1=m1 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.1$marginals[[1]])
m2.clone = inla.expectation(function(x) 1/x, seeds.hyperpar.clone.1$marginals[[1]])
stdev.clone.1 = sqrt(k*(m2.clone- mean.clone.1^2))

####---------------------------------------------------------------------------
#### Fitting by HDC method with flat prior ######

formula.clone.1.2 = r.k ~ x1.k+x2.k+f(plate.k,model="iid")
mod.seeds.clone.1.2 = inla(formula.clone.1.2,data=clone.data,
control.family=list(prior="flat"), family="binomial",Ntrials=n.k,
control.fixed=list(mean=c(1,1),prec=c(0.001,0.001)))

seeds.hyperpar.clone.1.2 = inla.hyperpar(mod.seeds.clone.1.2)

#### Fixed effects summaries
summary(mod.seeds.clone.1.2)
sd.clone.1.2<-as.vector(sqrt(k)*mod.seeds.clone.1.2$summary.fixed[,2])

#### Preparing elements for plotting 
marginal.clone.2.intercept<-mod.seeds.clone.1.2$marginals.fixed[[1]]
marginal.clone.2.x1<-mod.seeds.clone.1.2$marginals.fixed$x1
marginal.clone.2.x2<-mod.seeds.clone.1.2$marginals.fixed$x2
hyp.k.2<-seeds.hyperpar.clone.1.2$marginals[[1]]

#### Variance component summaries
mean.clone.1.2 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.1.2$marginals[[1]])
m2.clone.2 = inla.expectation(function(x) 1/x, seeds.hyperpar.clone.1.2$marginals[[1]])
stdev.clone.1.2 = sqrt(k*(m2.clone.2- mean.clone.1.2^2))

####---------------------------------------------------------------------------
#### Fitting by HDC method with vague (gamma) prior ######

formula.clone.1.3 = r.k ~ x1.k+x2.k+f(plate.k,model="iid",param=c(.01,0.01))
mod.seeds.clone.1.3 = inla(formula.clone.1.3,data=clone.data,
family="binomial",Ntrials=n.k, control.fixed=list(mean=c(-2,-1),prec=c(0.1,0.1)))
seeds.hyperpar.clone.1.3 = inla.hyperpar(mod.seeds.clone.1.3)

#### Fixed effects summaries
summary(mod.seeds.clone.1.3)
sd.clone.1.3<-as.vector(sqrt(k)*mod.seeds.clone.1.3$summary.fixed[,2])
inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.1.3$marginals[[1]])

#### Preparing elements for plotting 
marginal.clone.3.intercept<-mod.seeds.clone.1.3$marginals.fixed[[1]]
marginal.clone.3.x1<-mod.seeds.clone.1.3$marginals.fixed$x1
marginal.clone.3.x2<-mod.seeds.clone.1.3$marginals.fixed$x2
hyp.k.3<-seeds.hyperpar.clone.1.3$marginals[[1]]

#### Variance component summaries
mean.clone.1.3=m1.3 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.1.3$marginals[[1]])
m2.clone.3 = inla.expectation(function(x) 1/x, seeds.hyperpar.clone.1.3$marginals[[1]])
stdev.clone.1.3 = sqrt(k*(m2.clone.3- mean.clone.1.3^2))

####---------------------------------------------------------------------------
###################################################
#### Plotting the results to compare HDC-based densities with respect to priors

par(mfrow=c(2,2))

plot(marginal.clone.intercept,main=expression(paste("HDC-based density of",~ beta[0])),
col=1,xlab="",ylab="",lty=1)
lines(marginal.clone.2.intercept,lty=3,col=3)
lines(marginal.clone.3.intercept,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.clone.x1,main=expression(paste("DC-based density of",~ beta[1])),
col=1,xlab="",ylab="",lty=1)
lines(marginal.clone.2.x1,lty=3,col=3)
lines(marginal.clone.3.x1,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.clone.x2,main=expression(paste("DC-based density of",~ beta[2])),
col=1,xlab="",ylab="",lty=1)
lines(marginal.clone.2.x2,lty=3,col=3)
lines(marginal.clone.3.x2,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(hyp.k.1,main=expression(paste("DC-based density of",~ sigma^{-2})),
col=1,xlab="",ylab="",lty=1,ylim=c(0,.7))
lines(hyp.k.2,lty=3,col=3)
lines(hyp.k.3,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

###############################################################################
############### Second Model for Interaction Model ###########################

####----------------------------------------------------------------
### Fitting by AGHQ method with 15 nodes #####

(interaction.effects <- glmer(r/n ~ x1 + x2 + I(x1*x2) + (1 | plate) , 
			family = "binomial", data = Seeds, nAGQ=15))

####----------------------------------------------------------------
## Fitting by INLA
seeds.inla.fit.2 = inla(r ~ x1 + x2+I(x1*x2) + f(plate, model="iid",
        param=c(.5, .0164)), data=Seeds, family="binomial", Ntrials=n )

seeds.hyperpar.2 = inla.hyperpar(seeds.inla.fit.2)

### Fixed effects summaries
summary(seeds.inla.fit.2)

### Variance components summaries
mean.2 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.2$marginals[[1]])
m2 = inla.expectation(function(x) 1/x, seeds.hyperpar.2$marginals[[1]])
stdev.2 = sqrt(m2- mean.2^2)

####---------------------------------------------------------------------------
################ Fitting by HDC with informative prior

formula.clone.int.1 = r.k ~ x1.k+x2.k+I(x1.k*x2.k)+f(plate.k,model="iid",param=c(.5,.0164))
mod.seeds.clone.int.1 = inla(formula.clone.int.1,data=clone.data,family="binomial",Ntrials=n.k,
control.fixed=list(mean=c(1,1,1),prec=c(0.001,0.001,0.001)))

seeds.hyperpar.clone.int.1 = inla.hyperpar(mod.seeds.clone.int.1)

### Fixed effects summaries
summary(mod.seeds.clone.int.1)
sd.clone.int.1<-as.vector(sqrt(k)*mod.seeds.clone.int.1$summary.fixed[,2])

#### Preparing elements for plotting 
marginal.clone.int.1.intercept<-mod.seeds.clone.int.1$marginals.fixed[[1]]
marginal.clone.int.1.x1<-mod.seeds.clone.int.1$marginals.fixed$x1.k
marginal.clone.int.1.x2<-mod.seeds.clone.int.1$marginals.fixed$x2.k
marginal.clone.int.1.interaction<-mod.seeds.clone.int.1$marginals.fixed[[4]]
hyp.k.int.1<-seeds.hyperpar.clone.int.1$marginals$`Precision for plate.k`

### Variance components summaries
mean.clone.int.1 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.int.1$marginals[[1]])
m2.clone.int.1 = inla.expectation(function(x) 1/x, seeds.hyperpar.clone.int.1$marginals[[1]])
stdev.clone.int.1 = sqrt(k*(m2.clone.int.1- mean.clone.int.1^2))


# -----------------------------------------------------------
################ Fitting by HDC with flat prior

formula.clone.int.2 = r.k ~ x1.k+x2.k+I(x1.k*x2.k)+f(plate.k,model="iid")
mod.seeds.clone.int.2 = inla(formula.clone.int.2,data=clone.data,
control.family=list(prior="flat"), family="binomial",Ntrials=n.k,
control.fixed=list(mean=c(-1,-2,-1),prec=c(0.1,0.1,0.1)))

seeds.hyperpar.clone.int.2 = inla.hyperpar(mod.seeds.clone.int.2)

### Fixed effects summaries
summary(mod.seeds.clone.int.2)
sd.clone.int.2<-as.vector(sqrt(k)*mod.seeds.clone.int.2$summary.fixed[,2])

#### Preparing elements for plotting 
marginal.clone.int.2.intercept<-mod.seeds.clone.int.2$marginals.fixed[[1]]
marginal.clone.int.2.x1<-mod.seeds.clone.int.2$marginals.fixed$x1.k
marginal.clone.int.2.x2<-mod.seeds.clone.int.2$marginals.fixed$x2.k
marginal.clone.int.2.interaction<-mod.seeds.clone.int.2$marginals.fixed[[4]]
hyp.k.int.2<-seeds.hyperpar.clone.int.2$marginals$`Precision for plate.k`

### Variance components summaries
mean.clone.int.2 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.int.2$marginals[[1]])
m2.clone.int.2 = inla.expectation(function(x) 1/x, seeds.hyperpar.clone.int.2$marginals[[1]])
stdev.clone.int.2 = sqrt(k*(m2.clone.int.2- mean.clone.int.2^2))

# -----------------------------------------------------------
################ Fitting by HDC with vague (gamma) prior

formula.clone.int.3 = r.k ~ x1.k+x2.k+I(x1.k*x2.k)+f(plate.k,model="iid",param=c(.01,0.01))
mod.seeds.clone.int.3 = inla(formula.clone.int.3,data=clone.data,
family="binomial",Ntrials=n.k, control.fixed=list(mean=c(-2,2,1),prec=c(0.001,0.001,0.001)))

seeds.hyperpar.clone.int.3 = inla.hyperpar(mod.seeds.clone.int.3)

### Fixed effects summaries
summary(mod.seeds.clone.int.3)
sd.clone.int.3<-as.vector(sqrt(k)*mod.seeds.clone.int.3$summary.fixed[,2])
inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.int.3$marginals[[1]])

#### Preparing elements for plotting
marginal.clone.int.3.intercept<-mod.seeds.clone.int.3$marginals.fixed[[1]]
marginal.clone.int.3.x1<-mod.seeds.clone.int.3$marginals.fixed$x1.k
marginal.clone.int.3.x2<-mod.seeds.clone.int.3$marginals.fixed$x2.k
marginal.clone.int.3.interaction<-mod.seeds.clone.int.3$marginals.fixed[[4]]
hyp.k.int.3<-seeds.hyperpar.clone.int.3$marginals$`Precision for plate.k`

## Variance components summaries
mean.clone.int.3 = inla.expectation(function(x) 1/x^.5, seeds.hyperpar.clone.int.3$marginals[[1]])
m2.clone.int.3 = inla.expectation(function(x) 1/x, seeds.hyperpar.clone.int.3$marginals[[1]])
stdev.clone.int.3 = sqrt(k*(m2.clone.int.3- mean.clone.int.3^2))

# -----------------------------------------------------------
#### Plotting the results 

par(mfrow=c(2,2))

plot(marginal.clone.int.1.intercept,
main=expression(paste("DC-based distribution",~ beta[0])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.clone.int.2.x1,lty=3,col=3)
lines(marginal.clone.int.3.x1,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.clone.int.1.x1,
main=expression(paste("DC-based distribution",~ beta[1])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.clone.int.2.x1,lty=3,col=3)
lines(marginal.clone.int.3.x1,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.clone.int.1.x2,
main=expression(paste("DC-based distribution",~ beta[2])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.clone.int.2.x2,lty=3,col=3)
lines(marginal.clone.int.3.x2,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.clone.int.1.interaction,
main=expression(paste("DC-based distribution",~ beta[3])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.clone.int.2.interaction,lty=3,col=3)
lines(marginal.clone.int.3.interaction,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(hyp.k.int.1,main=expression(paste("DC-based distribution",~ sigma^{-2})),
col=1,xlab="",ylab="Density",lty=1,type="l",ylim=c(0,.4))
lines(hyp.k.int.2,lty=3,col=3)
lines(hyp.k.int.3,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

# -------------------------------------------------------------------------
#### Comparing DC and HDC results and their computing times ####

formula = r ~ x1+x2+f(plate,model="iid",param=c(.5,.0164))
formula.clone = r.k ~ x1.k+x2.k+f(plate.k,model="iid",param=c(.5,.0164))

## Run HDC
result.inla = inla(formula.clone, data=clone.data, family="binomial",
        Ntrials=n.k, verbose=T)

## Run DC

wd.dc = tempfile()
## same as before
result.dc = system.time(inla(formula.clone,data=clone.data,
	family="binomial",Ntrials=n.k,
	working.directory = wd.dc,
  	keep = TRUE,
 	inla.arg = "-m mcmc -N 10000 -T 10 -S 1"))

## Compare cpu times
time.inla = result.inla$cpu.used
time.inla
result.dc

## compare curves  (for example densities for the hyperparameters)
hyp.inla = result.inla$marginals.hyperpar[[1]]
hyp.dc = scan(paste(wd.dc,"/results.files/hyperparameter-random.effect00000001-parameter-user-scale/trace.dat",sep=""))
plot(hyp.inla, type="l",lty=1,col=1,xlab="",ylab="",
main=expression(paste("hybrid DC and DC-based densities of"~ sigma^{-2})),ylim=c(0,0.5))
lines(density(hyp.dc),lty=2,col=2)

## compare curve for fixed effects (for example \beta_1)
int.inla = result.inla$marginals.fixed[[2]]
int.dc = scan(paste(wd.dc,"/results.files/fixed.effect00000002/trace.dat",sep=""))
plot(int.inla,type="l",lty=1,col=1,ylim=c(0,20),xlab="",ylab="",
main=expression(paste("hybrid DC and DC based densities of",~ beta[1])))
lines(density(int.dc),col=2,lty=2)

############################################################################
# -------------------------------------------------------------------------
#########################################################################
###############     Longitudianl Data Example: Epilepsy Data    #########

### Loading data 
## Visit is created by glmmAK, it corresponds to Breslow and Clayton's
## Visit/10, because the codes are -3,-1,1,3.
require (glmmAK)
data(epilepticBC)
epil = epilepticBC
epil$id2=epil$id
epil$rand=1:nrow(epil)
epil$V4=epil$visit==4
epil$newid=rep(1:(nrow(epil)/4), each=4)

###---------------------------------------------------------------------
### Constructing cloned data set
####################################

epil.data<-epil
patient=59
n1=patient*4
k=100
n=n1*k

"id.k"<-rep(epil.data[,1],k)
"visit.k"<-rep(epil.data[,2],k)
"seizure0.k"<-rep(epil.data[,3],k)
"age.k"<-rep(epil.data[,4],k)
"Seizure.k"<-rep(epil.data[,5],k)
"Base.k"<-rep(epil.data[,6],k)
"Trt.k"<-rep(epil.data[,7],k)
"Base.Trt.k"<-rep(epil.data[,8],k)
"Age.k"<-rep(epil.data[,9],k)
"Visit.k"<-rep(epil.data[,10],k)
"id2.k"<-rep(epil.data[,11],k)
"rand.k"<-rep(epil.data[,12],k)
"V4.k"<-rep(epil.data[,13],k)
"newid.k"<-rep(epil.data[,14],k)


clone.data<-cbind(id.k,visit.k,seizure0.k,age.k,Seizure.k,
Base.k,Trt.k,Base.Trt.k,Age.k,Visit.k,id2.k,rand.k,V4.k,newid.k)

clone.data[,12]<-1:n
clone.data[,1]<-clone.data[,11]<-rep((1+100):(patient*k+100),each=4)
clone.data[,2]<-clone.data[,14]<-rep(1:patient*k,each=4)
clone.data<-as.data.frame(clone.data)

##############################################################################
##################   First random intercept model for Epilepsy data  #########
##################   proposed by Fong et al. (2009)  #########################
###---------------------------------------------------------------------

### Fitting by AGHQ method with 15 nodes 

(first.model <- glmer(Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 +
(1 | id) , family = "poisson", data = epil, nAGQ=15))

###---------------------------------------------------------------------

### Fitting by INLA method 

formula=Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 +
        f(id,model="iid",param=c(2, 1.140),diagonal=0)

epil.inla.fit.1 = inla(formula, data=epil, family="poisson" ,
			control.compute=list(hyperpar=T,dic=T))

epil.hyperpar.1 = inla.hyperpar(epil.inla.fit.1)

### Fixed effects summaries
summary(epil.inla.fit.1)

### Variance component summaries
mean.1=m1 = inla.expectation(function(x) 1/x^.5, epil.hyperpar.1$marginals[[1]])
m2 = inla.expectation(function(x) 1/x, epil.hyperpar.1$marginals[[1]])
stdev.1 = sqrt(m2- mean.1^2)

###---------------------------------------------------------------------

### Fitting by HDC method with informative prior

epil.dc.fit.1 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k)
+Age.k + V4.k+f(id.k,model="iid",param=c(2, 1.140)), data=clone.data,
family="poisson" )

epil.dc.hyperpar.1 = inla.hyperpar(epil.dc.fit.1)

### Fixed effects summaries
summary(epil.dc.fit.1)
sd.clone.1<-as.vector(sqrt(k)*epil.dc.fit.1$summary.fixed[,2])

### Variance component summaries
mean.dc.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.1$marginals[[1]])
m2.dc.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.1$marginals[[1]])
stdev.clone.1 = sqrt(k*(m2.dc.1- mean.dc.1^2))

#### Preparing elements to plot the results
marginal.Intercept.1<-epil.dc.fit.1$marginals.fixed[[1]]
marginal.Base.1<-epil.dc.fit.1$marginals.fixed$Base.k
marginal.Trt.1<-epil.dc.fit.1$marginals.fixed$Trt.k
marginal.BaseTrt.1<-epil.dc.fit.1$marginals.fixed[[4]]
marginal.Age.1<-epil.dc.fit.1$marginals.fixed$Age.k
marginal.V4.1<-epil.dc.fit.1$marginals.fixed$V4.k
hyp.clone.1<-epil.dc.hyperpar.1$marginals$`Precision for id.k`

###---------------------------------------------------------------------

### Fitting by HDC method with flat prior

epil.dc.fit.2 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k)
+Age.k + V4.k+f(id.k,model="iid"), data=clone.data,control.family=list(prior="flat"),
family="poisson",control.fixed=list(mean=c(1,1,1,1,0),prec=c(0.001,0.001,0.01,0.001,
0.01)))

epil.dc.hyperpar.2 = inla.hyperpar(epil.dc.fit.2)

### Fixed effects summaries
summary(epil.dc.fit.2)
sd.clone.2<-as.vector(sqrt(k)*epil.dc.fit.2$summary.fixed[,2])

### Variance component summaries
mean.dc.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2$marginals[[1]])
m2.dc.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2$marginals[[1]])
stdev.clone.2 = sqrt(k*(m2.dc.2- mean.dc.2^2))

#### Preparing elements to plot the results
marginal.Intercept.2<-epil.dc.fit.2$marginals.fixed[[1]]
marginal.Base.2<-epil.dc.fit.2$marginals.fixed$Base.k
marginal.Trt.2<-epil.dc.fit.2$marginals.fixed$Trt.k
marginal.BaseTrt.2<-epil.dc.fit.2$marginals.fixed[[4]]
marginal.Age.2<-epil.dc.fit.2$marginals.fixed$Age.k
marginal.V4.2<-epil.dc.fit.2$marginals.fixed$V4.k
hyp.clone.2<-epil.dc.hyperpar.2$marginals$`Precision for id.k`


###---------------------------------------------------------------------

### Fitting by HDC method with vague (gamma) prior

epil.dc.fit.3 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k)
+Age.k + V4.k+f(id.k,model="iid",param=c(0.05,0.02)), data=clone.data,
family="poisson",control.fixed=list(mean=c(-1,-2,2,0,1),prec=c(0.01,0.001,0.001,0.01,
0.001)))

epil.dc.hyperpar.3 = inla.hyperpar(epil.dc.fit.3)

### Fixed effects summaries
summary(epil.dc.fit.3)
sd.clone.3<-as.vector(sqrt(k)*epil.dc.fit.3$summary.fixed[,2])

### Variance component summaries
mean.dc.3 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3$marginals[[1]])
m2.dc.3 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3$marginals[[1]])
stdev.clone.3 = sqrt(k*(m2.dc.3- mean.dc.3^2))

#### Preparing elements to plot the results
marginal.Intercept.3<-epil.dc.fit.3$marginals.fixed[[1]]
marginal.Base.3<-epil.dc.fit.3$marginals.fixed$Base.k
marginal.Trt.3<-epil.dc.fit.3$marginals.fixed$Trt.k
marginal.BaseTrt.3<-epil.dc.fit.3$marginals.fixed[[4]]
marginal.Age.3<-epil.dc.fit.3$marginals.fixed$Age.k
marginal.V4.3<-epil.dc.fit.3$marginals.fixed$V4.k
hyp.clone.3<-epil.dc.hyperpar.3$marginals$`Precision for id.k`

###---------------------------------------------------------------------
#### Plotting the results

par(mfrow=c(3,2))

plot(marginal.Intercept.1,main=expression(paste(beta[0])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Intercept.2,lty=3,col=3)
lines(marginal.Intercept.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))


plot(marginal.Base.1,main=expression(paste(beta[1])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Base.2,lty=3,col=3)
lines(marginal.Base.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Trt.1,main=expression(paste(beta[2])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Trt.2,lty=3,col=3)
lines(marginal.Trt.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.BaseTrt.1,main=expression(paste(beta[3])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.BaseTrt.2,lty=3,col=3)
lines(marginal.BaseTrt.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Age.1,main=expression(paste(beta[4])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Age.2,lty=3,col=3)
lines(marginal.Age.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.V4.1,main=expression(paste(beta[5])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.V4.2,lty=3,col=3)
lines(marginal.V4.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(hyp.clone.1,main=expression(paste("HDC-based distribution of",~ sigma^{-2})),
col=1,xlab="",ylab="",lty=1,type="l")
lines(hyp.clone.2,lty=3,col=3)
lines(hyp.clone.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))


##############################################################################
##################   Second random intercept with measurement error  #########
############     model for Epilepsy data  proposed by Fong et al. (2009)  ####
###---------------------------------------------------------------------

### Fitting by AGHQ method with 15 nodes 

(second.model<- glmer(Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 +
(1 | id) + (1|rand), family = "poisson", data = epil))
###---------------------------------------------------------------------

### Fitting by INLA method 

formula.2=Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 +
        f(id,model="iid",param=c(2, 1.240)) +
        f(rand,model="iid",param=c(2, 1.140), diagonal=0)

epil.inla.fit.2 = inla(formula.2, data=epil,
        family="poisson" )

epil.hyperpar.2 = inla.hyperpar(epil.inla.fit.2)

### Fixed effects summaries
summary(epil.inla.fit.2)

### Variance component summaries
mean.1 = inla.expectation(function(x) 1/x^.5, epil.hyperpar.2$marginals[[1]])
m1 = inla.expectation(function(x) 1/x, epil.hyperpar.2$marginals[[1]])
stdev.1 = sqrt(m1- mean.1^2)
mean.2=m1 = inla.expectation(function(x) 1/x^.5, epil.hyperpar.2$marginals[[2]])
m2 = inla.expectation(function(x) 1/x, epil.hyperpar.2$marginals[[2]])
stdev.2 = sqrt(m2- mean.2^2)


###---------------------------------------------------------------------

### Fitting by HDC method with informative prior

formula.clone.2=Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k) + Age.k + V4.k +
      f(id.k,model="iid",param=c(2, 1.240))+
	f(rand.k,model="iid",param=c(2, 1.140))
epil.dc.fit.2.1 = inla(formula.clone.2,data=clone.data,family="poisson" )

epil.dc.hyperpar.2.1 = inla.hyperpar(epil.dc.fit.2.1)

### Fixed effects summaries
summary(epil.dc.fit.2.1)
sd.clone.2.1<-as.vector(sqrt(k)*epil.dc.fit.2.1$summary.fixed[,2])

### Variance component summaries
mean.dc.2.1.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2.1$marginals[[1]])
m2.dc.2.1.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2.1$marginals[[1]])
stdev.clone.2.1.1 = sqrt(k*(m2.dc.2.1.1- mean.dc.2.1.1^2))

mean.dc.2.1.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2.1$marginals[[2]])
m2.dc.2.1.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2.1$marginals[[2]])
stdev.clone.2.1.2 = sqrt(k*(m2.dc.2.1.2- mean.dc.2.1.2^2))

#### Preparing elements to plot the results
marginal.Intercept.2.1<-epil.dc.fit.2.1$marginals.fixed[[1]]
marginal.Base.2.1<-epil.dc.fit.2.1$marginals.fixed$Base.k
marginal.Trt.2.1<-epil.dc.fit.2.1$marginals.fixed$Trt.k
marginal.BaseTrt.2.1<-epil.dc.fit.2.1$marginals.fixed[[4]]
marginal.Age.2.1<-epil.dc.fit.2.1$marginals.fixed$Age.k
marginal.V4.2.1<-epil.dc.fit.2.1$marginals.fixed$V4.k
hyp.clone.2.1.1<-epil.dc.hyperpar.2.1$marginals$`Precision for id.k`
hyp.clone.2.1.2<-epil.dc.hyperpar.2.1$marginals$`Precision for rand.k`

###---------------------------------------------------------------------

### Fitting by HDC method with flat prior

epil.dc.fit.2.2 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k) + Age.k + V4.k +
        f(id.k,model="iid")+f(rand.k,model="iid"),control.family=list(prior="flat"),
	  data=clone.data,family="poisson",
control.fixed=list(mean=c(1,1,1,1,0),prec=c(0.001,0.001,0.01,0.001,0.01)))

epil.dc.hyperpar.2.2 = inla.hyperpar(epil.dc.fit.2.2)

### Fixed effects summaries
summary(epil.dc.fit.2.2)
sd.clone.2.2<-as.vector(sqrt(k)*epil.dc.fit.2.2$summary.fixed[,2])

### Variance component summaries
mean.dc.2.2.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2.2$marginals[[1]])
m2.dc.2.2.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2.2$marginals[[1]])
stdev.clone.2.2.1 = sqrt(k*(m2.dc.2.2.1- mean.dc.2.2.1^2))

mean.dc.2.2.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2.2$marginals[[2]])
m2.dc.2.2.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2.2$marginals[[2]])
stdev.clone.2.2.2 = sqrt(k*(m2.dc.2.2.2- mean.dc.2.2.2^2))


#### Preparing elements to plot the results
marginal.Intercept.2.2<-epil.dc.fit.2.2$marginals.fixed[[1]]
marginal.Base.2.2<-epil.dc.fit.2.2$marginals.fixed$Base.k
marginal.Trt.2.2<-epil.dc.fit.2.2$marginals.fixed$Trt.k
marginal.BaseTrt.2.2<-epil.dc.fit.2.2$marginals.fixed[[4]]
marginal.Age.2.2<-epil.dc.fit.2.2$marginals.fixed$Age.k
marginal.V4.2.2<-epil.dc.fit.2.2$marginals.fixed$V4.k
hyp.clone.2.2.1<-epil.dc.hyperpar.2.2$marginals$`Precision for id.k`
hyp.clone.2.2.2<-epil.dc.hyperpar.2.2$marginals$`Precision for rand.k`


###---------------------------------------------------------------------

### Fitting by HDC method with vague (gamma) prior

epil.dc.fit.2.3 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k) + Age.k + V4.k +
        f(id.k,model="iid",param=c(0.05, 0.02))+f(rand.k,model="iid",param=c(0.01, 0.01)),
	  data=clone.data,family="poisson",
control.fixed=list(mean=c(-1,-2,2,0,1),prec=c(0.01,0.001,0.001,0.01,0.001)))

epil.dc.hyperpar.2.3 = inla.hyperpar(epil.dc.fit.2.3)

### Fixed effects summaries
summary(epil.dc.fit.2.3)
sd.clone.2.3<-as.vector(sqrt(k)*epil.dc.fit.2.3$summary.fixed[,2])

### Variance component summaries
mean.dc.2.3.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2.3$marginals[[1]])
m2.dc.2.3.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2.3$marginals[[1]])
stdev.clone.2.3.1 = sqrt(k*(m2.dc.2.3.1- mean.dc.2.3.1^2))

mean.dc.2.3.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.2.3$marginals[[2]])
m2.dc.2.3.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.2.3$marginals[[2]])
stdev.clone.2.3.2 = sqrt(k*(m2.dc.2.3.2- mean.dc.2.3.2^2))

#### Preparing elements to plot the results
marginal.Intercept.2.3<-epil.dc.fit.2.3$marginals.fixed[[1]]
marginal.Base.2.3<-epil.dc.fit.2.3$marginals.fixed$Base.k
marginal.Trt.2.3<-epil.dc.fit.2.3$marginals.fixed$Trt.k
marginal.BaseTrt.2.3<-epil.dc.fit.2.3$marginals.fixed[[4]]
marginal.Age.2.3<-epil.dc.fit.2.3$marginals.fixed$Age.k
marginal.V4.2.3<-epil.dc.fit.2.3$marginals.fixed$V4.k
hyp.clone.2.3.1<-epil.dc.hyperpar.2.3$marginals$`Precision for id.k`
hyp.clone.2.3.2<-epil.dc.hyperpar.2.3$marginals$`Precision for rand.k`

###---------------------------------------------------------------------
#### Plotting the results

par(mfrow=c(3,2))

plot(marginal.Intercept.2.1,main=expression(paste(beta[0])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Intercept.2.2,lty=3,col=3)
lines(marginal.Intercept.2.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Base.2.1,main=expression(paste(beta[1])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Base.2.2,lty=3,col=3)
lines(marginal.Base.2.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Trt.2.1,main=expression(paste("beta[2])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Trt.2.2,lty=3,col=3)
lines(marginal.Trt.2.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.BaseTrt.2.1,main=expression(paste(beta[3])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.BaseTrt.2.2,lty=3,col=3)
lines(marginal.BaseTrt.2.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Age.2.1,main=expression(paste(beta[4])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Age.2.2,lty=3,col=3)
lines(marginal.Age.2.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.V4.2.1,main=expression(paste(beta[5])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.V4.2.2,lty=3,col=3)
lines(marginal.V4.2.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

par(mfrow=c(2,1))
plot(hyp.clone.2.1.1,main=expression(paste("HDC-based distribution of",~ sigma[1]^{-2})),
col=1,xlab="",ylab="",lty=1,type="l")
lines(hyp.clone.2.2.1,lty=3,col=3)
lines(hyp.clone.2.3.1,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(hyp.clone.2.1.2,main=expression(paste("HDC-based distribution of",~ sigma[2]^{-2})),
col=1,xlab="",ylab="",lty=1,type="l")
lines(hyp.clone.2.2.2,lty=3,col=3)
lines(hyp.clone.2.3.2,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

##############################################################################
##################   Second random slope model  ##############################
############     for Epilepsy data  proposed by Fong et al. (2009)  ##########
###---------------------------------------------------------------------

### Fitting by AGHQ method with 15 nodes 

(third.model<- glmer(Seizure ~ Base + Trt + I(Base*Trt) + Age + Visit +
(Visit|id2), family = "poisson", nAGQ=15, data = epil))

###---------------------------------------------------------------------

### Fitting by INLA method 

epil.inla.fit.3 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age +
        Visit +f(id, model="2diidwishartpart0", param=c(5, 2.277904,
        1.692047, 0), diagonal=0) + f(id2, Visit,
        model="2diidwishartpart1", diagonal = 0), data=epil,
        family="poisson" )

epil.hyperpar.3 = inla.hyperpar(epil.inla.fit.3)

### Fixed effects summaries
summary(epil.inla.fit.3)

### Variance component summaries
mean.3=m1 = inla.expectation(function(x) 1/x^.5, epil.hyperpar.3$marginals[[1]])
m2 = inla.expectation(function(x) 1/x, epil.hyperpar.3$marginals[[1]])
stdev.3.1 = sqrt(m2- mean.3^2)

###---------------------------------------------------------------------

### Fitting by HDC method with first prior set

epil.dc.fit.3.1 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k) + Age.k +
        Visit.k +f(id.k, model="2diidwishartpart0", param=c(5, 2.277904,
        1.692047, 0), diagonal=0) + f(id2.k, Visit.k,
        model="2diidwishartpart1", diagonal = 0), data=clone.data,
        family="poisson" )

epil.dc.hyperpar.3.1 = inla.hyperpar(epil.dc.fit.3.1)

### Fixed effects summaries
summary(epil.dc.fit.3.1)
sd.clone.3.1<-as.vector(sqrt(k)*epil.dc.fit.3.1$summary.fixed[,2])

### Variance component summaries
mean.dc.3.1.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3.1$marginals[[1]])
m2.dc.3.1.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3.1$marginals[[1]])
stdev.clone.3.1.1 = sqrt(k*(m2.dc.3.1.1- mean.dc.3.1.1^2))

mean.dc.3.1.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3.1$marginals[[2]])
m2.dc.3.1.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3.1$marginals[[2]])
stdev.clone.3.1.2 = sqrt(k*(m2.dc.3.1.2- mean.dc.3.1.2^2))

#### Preparing elements to plot the results
marginal.Intercept.3.1<-epil.dc.fit.3.1$marginals.fixed[[1]]
marginal.Base.3.1<-epil.dc.fit.3.1$marginals.fixed$Base.k
marginal.Trt.3.1<-epil.dc.fit.3.1$marginals.fixed$Trt.k
marginal.BaseTrt.3.1<-epil.dc.fit.3.1$marginals.fixed[[4]]
marginal.Age.3.1<-epil.dc.fit.3.1$marginals.fixed$Age.k
marginal.V4.3.1<-epil.dc.fit.3.1$marginals.fixed[[6]]
hyp.clone.3.1.1<-epil.dc.hyperpar.3.1$marginals[[1]]
hyp.clone.3.1.2<-epil.dc.hyperpar.3.1$marginals[[2]]

###---------------------------------------------------------------------

### Fitting by HDC method with second prior set

epil.dc.fit.3.2 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k) + Age.k +
        Visit.k +f(id.k, model="2diidwishartpart0", param=c(4, 3,
        4, 0), diagonal=0) + f(id2.k, Visit.k,
        model="2diidwishartpart1", diagonal = 0), data=clone.data,
        family="poisson" )

epil.dc.hyperpar.3.2 = inla.hyperpar(epil.dc.fit.3.2)

### Fixed effects summaries
summary(epil.dc.fit.3.2)
sd.clone.3.2<-as.vector(sqrt(k)*epil.dc.fit.3.2$summary.fixed[,2])

### Variance component summaries
mean.dc.3.2.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3.2$marginals[[1]])
m2.dc.3.2.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3.2$marginals[[1]])
stdev.clone.3.2.1 = sqrt(k*(m2.dc.3.2.1- mean.dc.3.2.1^2))

mean.dc.3.2.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3.2$marginals[[2]])
m2.dc.3.2.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3.2$marginals[[2]])
stdev.clone.3.2.2 = sqrt(k*(m2.dc.3.2.2- mean.dc.3.2.2^2))


#### Preparing elements to plot the results
marginal.Intercept.3.2<-epil.dc.fit.3.2$marginals.fixed[[1]]
marginal.Base.3.2<-epil.dc.fit.3.2$marginals.fixed$Base.k
marginal.Trt.3.2<-epil.dc.fit.3.2$marginals.fixed$Trt.k
marginal.BaseTrt.3.2<-epil.dc.fit.3.2$marginals.fixed[[4]]
marginal.Age.3.2<-epil.dc.fit.3.2$marginals.fixed$Age.k
marginal.V4.3.2<-epil.dc.fit.3.2$marginals.fixed[[6]]
hyp.clone.3.2.1<-epil.dc.hyperpar.3.2$marginals[[1]]
hyp.clone.3.2.2<-epil.dc.hyperpar.3.2$marginals[[2]]


###---------------------------------------------------------------------

### Fitting by HDC method with third prior set

epil.dc.fit.3.3 = inla(Seizure.k ~ Base.k + Trt.k + I(Base.k*Trt.k) + Age.k +
        Visit.k +f(id.k, model="2diidwishartpart0", param=c(6, 0.5,
        0.5, 0), diagonal=0) + f(id2.k, Visit.k,
        model="2diidwishartpart1", diagonal = 0), data=clone.data,
        family="poisson" )

epil.dc.hyperpar.3.3 = inla.hyperpar(epil.dc.fit.3.3)

### Fixed effects summaries
summary(epil.dc.fit.3.3)
sd.clone.3.3<-as.vector(sqrt(k)*epil.dc.fit.3.3$summary.fixed[,2]

### Variance component summaries
mean.dc.3.3.1 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3.3$marginals[[1]])
m2.dc.3.3.1 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3.3$marginals[[1]])
stdev.clone.3.3.1 = sqrt(k*(m2.dc.3.3.1- mean.dc.3.3.1^2))

mean.dc.3.3.2 = inla.expectation(function(x) 1/x^.5, epil.dc.hyperpar.3.3$marginals[[2]])
m2.dc.3.3.2 = inla.expectation(function(x) 1/x, epil.dc.hyperpar.3.3$marginals[[2]])
stdev.clone.3.3.2 = sqrt(k*(m2.dc.3.3.2- mean.dc.3.3.2^2))

#### Preparing elements to plot the results
marginal.Intercept.3.3<-epil.dc.fit.3.3$marginals.fixed[[1]]
marginal.Base.3.3<-epil.dc.fit.3.3$marginals.fixed$Base.k
marginal.Trt.3.3<-epil.dc.fit.3.3$marginals.fixed$Trt.k
marginal.BaseTrt.3.3<-epil.dc.fit.3.3$marginals.fixed[[4]]
marginal.Age.3.3<-epil.dc.fit.3.3$marginals.fixed$Age.k
marginal.V4.3.3<-epil.dc.fit.3.3$marginals.fixed[[6]]
hyp.clone.3.3.1<-epil.dc.hyperpar.3.3$marginals[[1]]
hyp.clone.3.3.2<-epil.dc.hyperpar.3.3$marginals[[2]]

###---------------------------------------------------------------------
#### Plotting the results

par(mfrow=c(3,2))

plot(marginal.Intercept.3.1,main=expression(paste(beta[0])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Intercept.3.2,lty=3,col=3)
lines(marginal.Intercept.3.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Base.3.1,main=expression(paste(beta[1])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Base.3.2,lty=3,col=3)
lines(marginal.Base.3.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Trt.3.1,main=expression(paste(beta[2])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Trt.3.2,lty=3,col=3)
lines(marginal.Trt.3.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.BaseTrt.3.1,main=expression(paste(beta[3])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.BaseTrt.3.2,lty=3,col=3)
lines(marginal.BaseTrt.3.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.Age.3.1,main=expression(paste(beta[4])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.Age.3.2,lty=3,col=3)
lines(marginal.Age.3.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.V4.3.1,main=expression(paste(beta[5])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.V4.3.2,lty=3,col=3)
lines(marginal.V4.3.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

par(mfrow=c(2,1))
plot(hyp.clone.3.1.1,main=expression(paste(sigma[1]^{-2})),
col=1,xlab="",ylab="",lty=1,type="l")
lines(hyp.clone.3.2.1,lty=3,col=3)
lines(hyp.clone.3.3.1,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(hyp.clone.3.1.2,main=expression(paste(sigma[2]^{-2})),
col=1,xlab="",ylab="",lty=1,type="l")
lines(hyp.clone.3.2.2,lty=3,col=3)
lines(hyp.clone.3.3.2,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))


############################################################################
# -------------------------------------------------------------------------
##################################################################################
########### Crossed Random Effects Data Example: Salamander Mating Data  #########

# -------------------------------------------------------------------------

### Loading data 

load("salam.RData")

## organize data into a form suitable for logistic regression
dat0=data.frame("y"=c(salam$y), "fW"=as.integer(salam$x[,"W/R"]==1 | salam$x[,"W/W"]==1), 
    "mW"=as.integer(salam$x[,"R/W"]==1 | salam$x[,"W/W"]==1), 
    "WW"=as.integer(salam$x[,"W/W"]==1 ) )
## add salamander id
id = t( apply(salam$z, 1, function(x) {
        tmp = which (x==1)
        tmp[2] = tmp[2] - 20
        tmp
    }) ) 
## ids are suitable for model A and C, but not B
id.modA = rbind(id, id+40, id+20)
colnames (id.modA) = c("f.modA","m.modA")
dat0=cbind (dat0, id.modA, group=1)
dat0$experiment=as.factor(rep(1:3, each=120))
dat0$group=as.factor(dat0$group)

salamander = dat0
salamander.e1 = subset (dat0, dat0$experiment==1)
salamander.e2 = subset (dat0, dat0$experiment==2)
salamander.e3 = subset (dat0, dat0$experiment==3)

### Constructing Cloned data set

k=100

dat0=salamander.e1
dat0$no<-1:120

dat0.1 = subset (dat0, dat0$no==(1:20))
dat0.2 = subset (dat0, dat0$no==(21:40))
dat0.3 = subset (dat0, dat0$no==(41:60))
dat0.4 = subset (dat0, dat0$no==(61:80))
dat0.5 = subset (dat0, dat0$no==(81:100))
dat0.6 = subset (dat0, dat0$no==(101:120))
###
data.clone<-matrix(NA,nrow=20*k,ncol=9)
cloning<-function(k){
for (i in 1:9){
	data.clone[,i]<-rep(dat0.1[,i],k)
	}
return(data.clone)
}
data.0.1.clone<-cloning(k)

data.0.1.clone[,5]<-1:(20*k)
for (i in 2:k){
data.0.1.clone[(((i-1)*20+1):(i*20)),6]<-data.0.1.clone[(((i-1)*20+1):(i*20)),6]+(i-1)*20
}
data.clone.1<-data.0.1.clone[,-9]
#####

data.clone<-matrix(NA,nrow=20*k,ncol=9)
cloning<-function(k){
for (i in 1:9){
	data.clone[,i]<-rep(dat0.2[,i],k)
	}
return(data.clone)
}
data.0.2.clone<-cloning(k)


data.0.2.clone[,5]<-1:(20*k)
for (i in 2:k){
data.0.2.clone[(((i-1)*20+1):(i*20)),6]<-data.0.2.clone[(((i-1)*20+1):(i*20)),6]+(i-1)*20
}
data.clone.2<-data.0.2.clone[,-9]
####

data.clone<-matrix(NA,nrow=20*k,ncol=9)
cloning<-function(k){
for (i in 1:9){
	data.clone[,i]<-rep(dat0.3[,i],k)
	}
return(data.clone)
}
data.0.3.clone<-cloning(k)


data.0.3.clone[,5]<-1:(20*k)
for (i in 2:k){
data.0.3.clone[(((i-1)*20+1):(i*20)),6]<-data.0.3.clone[(((i-1)*20+1):(i*20)),6]+(i-1)*20
}
data.clone.3<-data.0.3.clone[,-9]
#######

data.clone<-matrix(NA,nrow=20*k,ncol=9)
cloning<-function(k){
for (i in 1:9){
	data.clone[,i]<-rep(dat0.4[,i],k)
	}
return(data.clone)
}
data.0.4.clone<-cloning(k)


data.0.4.clone[,5]<-1:(20*k)
for (i in 2:k){
data.0.4.clone[(((i-1)*20+1):(i*20)),6]<-data.0.4.clone[(((i-1)*20+1):(i*20)),6]+(i-1)*20
}
data.clone.4<-data.0.4.clone[,-9]
####

data.clone<-matrix(NA,nrow=20*k,ncol=9)
cloning<-function(k){
for (i in 1:9){
	data.clone[,i]<-rep(dat0.5[,i],k)
	}
return(data.clone)
}
data.0.5.clone<-cloning(k)


data.0.5.clone[,5]<-1:(20*k)
for (i in 2:k){
data.0.5.clone[(((i-1)*20+1):(i*20)),6]<-data.0.5.clone[(((i-1)*20+1):(i*20)),6]+(i-1)*20
}
data.clone.5<-data.0.5.clone[,-9]
#####
data.clone<-matrix(NA,nrow=20*k,ncol=9)
cloning<-function(k){
for (i in 1:9){
	data.clone[,i]<-rep(dat0.6[,i],k)
	}
return(data.clone)
}
data.0.6.clone<-cloning(k)


data.0.6.clone[,5]<-1:(20*k)
for (i in 2:k){
data.0.6.clone[(((i-1)*20+1):(i*20)),6]<-data.0.6.clone[(((i-1)*20+1):(i*20)),6]+(i-1)*20
}
data.clone.6<-data.0.6.clone[,-9]
#############

data.clone<-rbind(data.clone.1,data.clone.2,data.clone.3,data.clone.4,
data.clone.5,data.clone.6)

aa<-data.clone<-as.data.frame(data.clone)

#########################################################################
###-----------------------------------------------------------------------
#########################   Summer Experiment ############################
#### For two other experiments the codes are the same ####################
###-----------------------------------------------------------------------
### Fitting by INLA approach

formula=y~fW+mW+WW + f(f.modA, model="iid", param=c(1,.622)) + 
f(m.modA, model="iid", param=c(1,.622))

salamander.e1.inla.fit = inla(formula, 
        family="binomial", data=salamander.e1, Ntrials=rep(1,nrow(salamander.e1)))
salamander.e1.hyperpar = inla.hyperpar (salamander.e1.inla.fit)

### Fixed effects summaries
summary(salamander.e1.inla.fit)

### Variance components summaries
mean.inla.1 = inla.expectation(function(x) 1/x^.5, salamander.e1.hyperpar$marginals[[1]])
mean.inla.2 = inla.expectation(function(x) 1/x^.5, salamander.e1.hyperpar$marginals[[2]])

m2.inla.1 = inla.expectation(function(x) 1/x, salamander.e1.hyperpar$marginals[[1]])
m2.inla.2 = inla.expectation(function(x) 1/x, salamander.e1.hyperpar$marginals[[2]])
stdev.inla.1 = sqrt(k*(m2.inla.1- mean.inla.1^2))
stdev.inla.2 = sqrt(k*(m2.inla.2- mean.inla.2^2))

###-----------------------------------------------------------------------
### Fitting by HDC approach with informative prior set 

formula.clone=aa[,1]~aa[,2]+aa[,3]+aa[,4]+f(aa[,5],model="iid",param=c(1,.622)) +
f(aa[,6], model="iid", param=c(1,.622))

salamander.e1.clone.1 = inla(formula.clone,
family="binomial", data=aa, Ntrials=rep(1,nrow(aa)))
salamander.clone.e1.hyperpar.1 = inla.hyperpar (salamander.e1.clone.1)

### Fixed effects summaries
summary(salamander.e1.clone.1)
sd.clone.e1<-as.vector(sqrt(k)*salamander.e1.clone.1$summary.fixed[,2])

### Variance components summaries
mean.clone.1.1=inla.expectation(function(x) 1/x^.5, salamander.clone.e1.hyperpar.1$marginals[[1]])
mean.clone.2.1=inla.expectation(function(x) 1/x^.5, salamander.clone.e1.hyperpar.1$marginals[[2]])

m2.clone.1.1 = inla.expectation(function(x) 1/x, salamander.clone.e1.hyperpar.1$marginals[[1]])
m2.clone.2.1 = inla.expectation(function(x) 1/x, salamander.clone.e1.hyperpar.1$marginals[[2]])
stdev.clone.1.1= sqrt(k*(m2.clone.1.1- mean.clone.1.1^2))
stdev.clone.2.1= sqrt(k*(m2.clone.2.1- mean.clone.2.1^2))

### Preparing elements to plot the results
marginal.intercept.1<-salamander.e1.clone.1$marginals.fixed[[1]]
marginal.F.WS.1<-salamander.e1.clone.1$marginals.fixed[[2]]
marginal.M.WS.1<-salamander.e1.clone.1$marginals.fixed[[3]]
marginal.FM.WS.1<-salamander.e1.clone.1$marginals.fixed[[4]]
hyp.F.1<-salamander.clone.e1.hyperpar.1$marginals[[1]]
hyp.M.1<-salamander.clone.e1.hyperpar.1$marginals[[2]]

###-----------------------------------------------------------------------
### Fitting by HDC approach with flat prior set 

salamander.e1.clone.2 = inla(aa[,1]~aa[,2]+aa[,3]+aa[,4]+f(aa[,5],model="iid")+ 
f(aa[,6], model="iid"), control.family=list(prior="flat"),
        family="binomial", data=aa, Ntrials=rep(1,nrow(aa)))
salamander.clone.e1.hyperpar.2 = inla.hyperpar (salamander.e1.clone.2)


### Fixed effects summaries
summary(salamander.e1.clone.2)
sd.clone.e2<-as.vector(sqrt(k)*salamander.e1.clone.2$summary.fixed[,2])

### Variance components summaries
mean.clone.1.2=inla.expectation(function(x) 1/x^.5, salamander.clone.e1.hyperpar.2$marginals[[1]])
mean.clone.2.2=inla.expectation(function(x) 1/x^.5, salamander.clone.e1.hyperpar.2$marginals[[2]])

m2.clone.1.2 = inla.expectation(function(x) 1/x, salamander.clone.e1.hyperpar.2$marginals[[1]])
m2.clone.2.2 = inla.expectation(function(x) 1/x, salamander.clone.e1.hyperpar.2$marginals[[2]])
stdev.clone.1.2= sqrt(k*(m2.clone.1.2- mean.clone.1.2^2))
stdev.clone.2.2= sqrt(k*(m2.clone.2.2- mean.clone.2.2^2))

### Preparing elements to plot the results
marginal.intercept.2<-salamander.e1.clone.2$marginals.fixed[[1]]
marginal.F.WS.2<-salamander.e1.clone.2$marginals.fixed[[2]]
marginal.M.WS.2<-salamander.e1.clone.2$marginals.fixed[[3]]
marginal.FM.WS.2<-salamander.e1.clone.2$marginals.fixed[[4]]
hyp.F.2<-salamander.clone.e1.hyperpar.2$marginals[[1]]
hyp.M.2<-salamander.clone.e1.hyperpar.2$marginals[[2]]

###-----------------------------------------------------------------------
### Fitting by HDC approach with vague (gamma) prior set 

salamander.e1.clone.3 = inla(aa[,1]~aa[,2]+aa[,3]+aa[,4]+f(aa[,5],model="iid",param=c(0.1,0.1)) 
+ f(aa[,6], model="iid", param=c(0.1,0.1)), 
        family="binomial", data=aa, Ntrials=rep(1,nrow(aa)))
salamander.clone.e1.hyperpar.3 = inla.hyperpar (salamander.e1.clone.3)

### Fixed effects summaries
summary(salamander.e1.clone.3)
sd.clone.e3<-as.vector(sqrt(k)*salamander.e1.clone.3$summary.fixed[,2])

### Variance components summaries
mean.clone.1.3=inla.expectation(function(x) 1/x^.5, salamander.clone.e1.hyperpar.3$marginals[[1]])
mean.clone.2.3=inla.expectation(function(x) 1/x^.5, salamander.clone.e1.hyperpar.3$marginals[[2]])

m2.clone.1.3 = inla.expectation(function(x) 1/x, salamander.clone.e1.hyperpar.3$marginals[[1]])
m2.clone.2.3 = inla.expectation(function(x) 1/x, salamander.clone.e1.hyperpar.3$marginals[[2]])
stdev.clone.1.3= sqrt(k*(m2.clone.1.3- mean.clone.1.3^2))
stdev.clone.2.3= sqrt(k*(m2.clone.2.3- mean.clone.2.3^2))


### Preparing elements to plot the results
marginal.intercept.3<-salamander.e1.clone.3$marginals.fixed[[1]]
marginal.F.WS.3<-salamander.e1.clone.3$marginals.fixed[[2]]
marginal.M.WS.3<-salamander.e1.clone.3$marginals.fixed[[3]]
marginal.FM.WS.3<-salamander.e1.clone.3$marginals.fixed[[4]]
hyp.F.3<-salamander.clone.e1.hyperpar.3$marginals[[1]]
hyp.M.3<-salamander.clone.e1.hyperpar.3$marginals[[2]]

###-----------------------------------------------------------------------

###############  Drawing HDC graphs ##################

par(mfrow=c(2,2))

plot(marginal.intercept.1,main=expression(paste("beta[0])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.intercept.2,lty=3,col=3)
lines(marginal.intercept.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.F.WS.1,main=expression(paste("beta[1])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.F.WS.2,lty=3,col=3)
lines(marginal.F.WS.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.M.WS.1,main=expression(paste(beta[2])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.M.WS.2,lty=3,col=3)
lines(marginal.M.WS.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(marginal.FM.WS.1,main=expression(paste(beta[3])),
col=1,xlab="",ylab="",lty=1,type="l")
lines(marginal.FM.WS.2,lty=3,col=3)
lines(marginal.FM.WS.3,lty=4,col=4)
#legend("topleft",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

par(mfrow=c(2,1))

plot(hyp.F.1,main=expression(paste(sigma[f]^{-2})),
col=1,xlab="",ylab="",lty=1,type="l")
lines(hyp.F.2,lty=3,col=3)
lines(hyp.F.3,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

plot(hyp.M.1,main=expression(paste(sigma[m]^{-2})),
col=1,xlab="",ylab="",lty=1,type="l",ylim=c(0,0.26))
lines(hyp.M.2,lty=3,col=3)
lines(hyp.M.3,lty=4,col=4)
#legend("topright",c("Priors 1","Priors 2","Priors 3")
#, bty="n",lty=c(1,3,4),col=c(1,3,4))

###----------------------------------------------------------------
################################################################################
########### Comparing computing times in DC and hybrid DC methods ################

wd.dc = tempfile()
## same as before
result.dc = system.time(inla(formula.clone,data=aa,
	family="binomial",Ntrials=rep(1,nrow(aa)),
	working.directory = wd.dc,
  	keep = TRUE,
 	inla.arg = "-m mcmc -N 10000 -T 10 -S 0.1"))

### computing cpu times
time.inla = salamander.e1.clone.1$cpu.used
time.inla
result.dc


