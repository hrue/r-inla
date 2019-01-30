##########################################################################################################
############################################################################################
#                                 Low Birth Weigth in Georgia (US)                         #                                  #
############################################################################################
#Import the data
data <- read.csv("data.final.csv")

#-- Prepare the map --#
library(maptools)
library(sp)
library(spdep)
georgia <- readShapePoly("co13_d00.shp")
## Need to drop extra polygons 98 (Macon county polygon), 100, 105 (Taylor county polygons) + 137 (Lee county polygon)
## These are very small and always adjacent to "main" polygon, so we can base neighborhood structure on "main" polygon:
rmIdx <- c(98, 100, 105, 137)
georgia <- georgia[-rmIdx,]
data.georgia = attr(georgia, "data")
#################################################
#Create the graph for adjacencies in INLA
#Need the non thinned sph file to do the adjacency matrix!!!
zzz <- poly2nb(georgia)
nb2INLA("Georgia.graph", zzz)
#this create a file called "LDN-INLA.adj" with the graph for INLA
Georgia.adj <- paste(getwd(),"/Georgia.graph",sep="")

#Order based on the map
order <- match(data.georgia$NAME,data[,1])
data<- data[order,]

#--Transform the data to be in the right format for INLA--#
low.vector <- as.vector(as.matrix(data[,2:12]))#by column
E.vector <- as.vector(as.matrix(data[,13:23]))#by column
year <- numeric(0)
for(i in 1:11){ 
  year<- append(year,rep(i,dim(data)[1]))
}
county<- as.factor(rep(data[,1],11))

data<- data.frame(y= low.vector, E= E.vector, ID.area=as.numeric(county), ID.area1=as.numeric(county), year=year,
                  ID.year = year, ID.year1=year, ID.area.year = seq(1,length(county)))
#############################################################
#--Prepare the model and run inla--#
library(INLA)
#Parametric model alpha + csii + (deltai + beta)*year
formula.ST1<- y ~ 1 + f(ID.area,model="bym",graph=Georgia.adj) +
  f(ID.area1,year,model="rw1") + (year-mean(year))
model.inla.ST1 <- inla(formula.ST1,family="poisson",data=data,E=E, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))

#Non Parametric model alpha + csii + gammaj + phij #No space time interaction yet!
#csii and are modelled through BYM
#gammaj are modelled as RW1
#phij are modelled as exchangeable
formula.ST2<- y ~ 1 + f(ID.area,model="bym",graph=Georgia.adj) +
  f(ID.year,model="rw1") + f(ID.year1,model="iid")
model.inla.ST2 <- inla(formula.ST2,family="poisson",data=data,E=E, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE))

#Non Parametric model alpha + csii + gammaj + phij + deltaij
#csii are modelled through BYM
#gammaj are modelled as RW1
#phij are modelled as exchangeable
#Interaction (deltaij) is modelled as exchangeable
formula.ST3<- y ~ 1 + f(ID.area,model="bym",graph=Georgia.adj) +
  f(ID.year,model="rw1") + f(ID.year1,model="iid") + f(ID.area.year,model="iid")

#To obtain the marginal of phij + gammaj we need to create the corresponding linear combinations and include these in the model 
lcs = inla.make.lincombs(ID.year = diag(11),  ID.year1 = diag(11))

model.inla.ST3 <- inla(formula.ST3,family="poisson",data=data,E=E, 
                       control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE),
                       lincomb=lcs,control.inla = list(lincomb.derived.only=TRUE))

#Put the temporal effect  (gammaj+phij) on the natural scale
temporal<-lapply(model.inla.ST3$marginals.lincomb.derived, function(X){
  marg <- inla.marginal.transform(function(x) exp(x), X)
  inla.emarginal(mean, marg)
})

#############################################################
#Computethe DIC as a tool for model choice
model.inla.ST1$dic$dic
model.inla.ST2$dic$dic
model.inla.ST3$dic$dic

#DIC components: Effective number of parameter (pd)
model.inla.ST1$dic$p.eff
model.inla.ST2$dic$p.eff
model.inla.ST3$dic$p.eff
#DIC components: mean.deviance
model.inla.ST1$dic$mean.deviance
model.inla.ST2$dic$mean.deviance
model.inla.ST3$dic$mean.deviance
######################################################
#The last model (with interaction) shows the best fit. Look at the results
#Obtain zetai exponentiating csii 
m <- model.inla.ST3$marginals.random[[1]][1:163]
zeta <- unlist(lapply(m,function(x)inla.emarginal(exp,x)))

#Probability that theta>1
a=0
inlaprob<-lapply(model.inla.ST3$marginals.random[[1]][1:163], function(X){
  1-inla.pmarginal(a, X)
})

Spatial.results<- data.frame(NAME=data.georgia$NAME,zeta=unlist(zeta), pp=unlist(inlaprob))

#Maps
#Create classes of SMRs
zeta.cutoff<- c(0.6, 0.9, 1.0, 1.1,1.8)
pp.cutoff <- c(0,0.2,0.8,1)
zeta=cut(Spatial.results$zeta,breaks=zeta.cutoff,include.lowest=TRUE)
pp=cut(Spatial.results$pp,breaks=pp.cutoff,include.lowest=TRUE)

maps.factors <- data.frame(NAME=data.georgia$NAME, zeta=zeta,pp=pp)
attr(georgia, "data")=data.frame(data.georgia, maps.factors)

trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia, zcol= "zeta", col.regions=gray(3.5:0.5/4),main="",par.settings=list(fontsize=list(text=17)))

trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia, zcol= "pp", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))

#Plot the National temporal trend
plot(seq(1,11),seq(0.8,1.2,length=11),type="n",xlab="year",ylab=expression(exp(gamma[t]+phi[t])))
lines(unlist(temporal))
abline(h=1,lty=2)
#######################
#Space-Time Interaction
delta <- data.frame(delta=model.inla.ST3$summary.random$ID.area.year[,2],year=data$ID.year,ID.area=data$ID.area)
delta.matrix <- matrix(delta[,1], 163,11,byrow=FALSE)
rownames(delta.matrix)<- delta[1:163,3]

#Space time probability>1
a=0
inlaprob.delta<-lapply(model.inla.ST3$marginals.random[[4]], function(X){
  1-inla.pmarginal(a, X)
})
pp.delta<-unlist(inlaprob.delta)

pp.cutoff.interaction <- c(0,0.2,0.8,1)
pp.delta.matrix <- matrix(pp.delta, 163,11,byrow=FALSE)
pp.delta.factor <- data.frame(NAME=data.georgia$NAME)
for(i in 1:11){
  pp.delta.factor.temp <- cut(pp.delta.matrix[,i],breaks=pp.cutoff.interaction,include.lowest=TRUE) 
  pp.delta.factor <- cbind(pp.delta.factor,pp.delta.factor.temp)
}
colnames(pp.delta.factor)<- c("NAME",seq(2000,2010))

#Maps
attr(georgia, "data")=data.frame(data.georgia, pp.delta.factor)

trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia, zcol="X2001", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia, zcol="X2004", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia, zcol="X2007", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia, zcol="X2010", col.regions=gray(2:0/2),main="",par.settings=list(fontsize=list(text=17)))
#############################################################
