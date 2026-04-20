############################################################################################                                  #
#                             London borough suicides                                      #
############################################################################################
## Edited 2023-09-14 to remove dependence on obsolete 'maptools' package                  ##
############################################################################################

#Load the package for building the map and import the shapefile
library(spdep)
london.gen = sf::st_read("../Data/LDNSuicides.shp")

#--- The data ---#
y=c(75,145,99,168,152,173,152,169,130,117,124,119,134,90,
    98,89,128,145,130,69,246,166,95,135,98,97,202,75,100,
    100,153,194)
E=c(80.7,169.8,123.2,139.5,169.1,107.2,179.8,160.4,147.5,
    116.8,102.8,91.8,119.6,114.8,131.1,136.1,116.6,98.5,
    88.8,79.8,144.9,134.7,98.9,118.6,130.6,96.1,127.1,97.7,
    88.5,121.4,156.8,114)
x1=c(0.87,-0.96,-0.84,0.13,-1.19,0.35,-0.84,-0.18,-0.39,0.74,
     1.93,0.24,0.59,-1.15,-0.8,-0.59,-0.12,1.43,-0.04,-1.24,1,
     0.53,-0.75,1.36,-0.93,-1.24,1.68,-1.04,2.38,0.03,-0.2,0.14)
x2=c(-1.02,-0.33,-1.43,-0.1,-0.98,1.77,-0.73,-0.2,-0.96,-0.58,
     0.36,1.48,0.46,-0.93,-1.62,-0.96,-0.48,0.81,2.41,-0.4,
     0.71,-0.05,-0.33,-0.47,-0.92,0.06,0.22,-0.73,0.1,-0.59,
     0.7,2.28)

names<- sort(london.gen$NAME)
data <- data.frame(NAME=names, y=y, E=E, x1=x1, x2=x2)
Nareas <- length(data[,1])

#In this section we create the adjacency graph
temp <- poly2nb(london.gen)
nb2INLA("LDN.graph", temp)

#This create a file called ``LDN-INLA.adj'' with the graph for INLA
LDN.adj <- paste(getwd(),"/LDN.graph",sep="")

#The order of the areas needs to be the same between the data and the spatial polygon object obtained importing the shapefile, so we re-order the data.
boroughs=london.gen
order <- match(boroughs$NAME,data$NAME)
data <- data[order,]
ID<-seq(1,32)
data <- cbind(ID,data)

#We now prepare the BYM model and run INLA
library(INLA)
formula <- y ~ 1 + f(ID, model="bym", graph=LDN.adj)
mod <- inla(formula,family="poisson",data=data,E=E)

#We calculate zeta=exp(csi) where csi=upsilon + nu
m <- mod$marginals.random$ID[1:Nareas]
zeta <- lapply(m,function(x)inla.emarginal(exp,x))

#We now calculate the probability that the spatial effects zeta are above 1,
#identifying areas with excess risk of suicides. This is equivalent to
#calculate the probability that csi is above 0, which is easier to obtain
a=0
inlaprob<-lapply(mod$marginals.random$ID[1:Nareas], function(X){
  1-inla.pmarginal(a, X)
})

#Finally we obtain the proportion of variance explained by the spatially structured component
#upsilon, taking the structured effect upsilon and calculating the empirical variance.
#First we create a matrix with rows equal to the number of areas and 1000 columns.
#Then for each area we extract 1000 values from the corresponding marginal distribution of upsilon
#and finally we calculate the empirical variance. We also extract the expected value of the variance
#for the unstructured component and build the spatial fractional variance.
mat.marg<-matrix(NA, nrow=Nareas, ncol=1000)
m<-mod$marginals.random$ID
for (i in 1:Nareas){
  u<-m[[Nareas+i]]
  s<-inla.rmarginal(1000, u)
  mat.marg[i,]<-s}
var.RRspatial<-mean(apply(mat.marg, 2, sd))^2
var.RRhet<-inla.emarginal(function(x) 1/x,
                          mod$marginals.hyper$"Precision for ID (iid component)")
var.RRspatial/(var.RRspatial+var.RRhet)
########################
#Now we add the covariates (deprivation - x1 and social fragmentation - x2) and repeat the steps
########################
formula.cov <- y ~ 1+ f(ID, model="bym", graph=LDN.adj) + x1 + x2
mod.cov <- inla(formula.cov,family="poisson",data=data,E=E)
mod.cov$summary.fixed
m <- mod.cov$marginals.random$ID[1:Nareas]
zeta.cov <- lapply(m,function(x)inla.emarginal(exp,x))
a=0
inlaprob.cov<-lapply(mod.cov$marginals.random$ID[1:Nareas], function(X){
  1-inla.pmarginal(a, X)
})
m<-mod.cov$marginals.random$ID
mat.marg<-matrix(NA, nrow=Nareas, ncol=1000)
for (i in 1:Nareas){
  u<-m[[Nareas+i]]
  s<-inla.rmarginal(1000, u)
  mat.marg[i,]<-s}
var.RRspatial<-mean(apply(mat.marg, 2, sd))^2
var.RRhet<-inla.emarginal(function(x) 1/x,
                          mod.cov$marginals.hyper$"Precision for ID (iid component)")
var.RRspatial/(var.RRspatial+var.RRhet)

#Finally we build the maps. First we create a dataset with all the relevant quantities
#and classes of SMR and posterior probabilities. Then transform the continuous SMR and
#posterior probabilities in factors, Merge the spatial polygon of London boroughs with the
#data and map the quantities.
Spatial.results<- data.frame(NAME=data$NAME,SMR=unlist(zeta),
                             pp=unlist(inlaprob), SMR.cov = unlist(zeta.cov), pp.cov = unlist(inlaprob.cov))
SMR.cutoff<- c(0.6, 0.9, 1.0, 1.1,  1.8)
pp.cutoff <- c(0,0.2,0.8,1)
#Transform SMR and pp in factors
SMR.DM=cut(Spatial.results$SMR,breaks=SMR.cutoff,include.lowest=TRUE)
pp.DM=cut(Spatial.results$pp,breaks=pp.cutoff,include.lowest=TRUE)
SMR.COV=cut(Spatial.results$SMR.cov,breaks=SMR.cutoff,include.lowest=TRUE)
pp.COV=cut(Spatial.results$pp.cov,breaks=pp.cutoff,include.lowest=TRUE)
maps.SMR.factors <- data.frame(NAME=data$NAME,
                               SMR.DM=SMR.DM,pp.DM=pp.DM,SMR.COV=SMR.COV,pp.COV=pp.COV)
boroughs <- merge(boroughs,maps.SMR.factors,by="NAME")

library(lattice)
boroughs <- sf::as_Spatial(boroughs)
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "SMR.DM", col.regions=gray(3.5:0.5/4),main="")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "pp.DM", col.regions=gray(2.5:0.5/3),main="")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "SMR.COV", col.regions=gray(3.5:0.5/4))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "pp.COV", col.regions=gray(2.5:0.5/3))
###########################################################################################
###########################################################################################
