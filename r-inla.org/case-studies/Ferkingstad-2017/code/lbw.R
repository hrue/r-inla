########################################################################################
## Code for "low birth weight" example, illustrating the use of the inla.cut()        ##
## function for group-wise model criticism; see ?inla.cut for documentation.          ## 
##                                                                                    ##
## For details see Section 4.2 of the following paper:                                ##
## Egil Ferkingstad, Leonhard Held and Havard Rue.                                    ##
## Fast and accurate Bayesian model criticism and conflict diagnostics using R-INLA.  ##
## Published in Stat, 6:331-344, 2017, doi: 10.1002/sta4.163.                         ##     
## Available as arXiv preprint arXiv:1708.03272: http://arxiv.org/abs/1708.03272      ##
########################################################################################

library(INLA)
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
#this creates a file called "Georgia.graph" with the graph for INLA
Georgia.adj <<- paste(getwd(),"/Georgia.graph",sep="")

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

data<- data.frame(y= low.vector, E= E.vector,
                  ID.area=as.numeric(county), ID.area1=as.numeric(county),
                  year=year,
                  ID.year = year, ID.year1=year,
                  ID.area.year = seq(1,length(county)))

data$intercept <- 1

#############################################################
#--Prepare the model and run inla--#

data$myear <- data$year - mean(data$year)
formula.ST1<- y ~ -1 + intercept + f(ID.area,model="bym",graph=Georgia.adj) + myear

rLBW <- inla(formula=formula.ST1,
                          data=data,
                          family="poisson",E=data$E,
                          control.predictor=list(link=1),
                          control.inla=list(h=0.05),
                     verbose=FALSE)

## Conflict detection with split in time:
pLBW.time <- inla.cut(rLBW, split.by="ID.year")

## plot results for split in time
## red line shows false discovery rate (FDR) limit for FDR=10%:
p <- pLBW.time
n.p <- length(p)
grid <- c(0:n.p)/n.p
y <- 0+0.1*grid
sel <- (sort(p) < y[-1])
par(pty="s")
plot(c(1:length(p))/n.p, sort(p), xlim=c(0,1), ylim=c(0,1), pch=19, cex=0.75,
     xlab="Rank of P-values divided by total number (30) of P-values", 
     ylab="Conflict P-values (ordered)", col=c(rep(2, sum(sel)), rep(1, n.p-sum(sel))))
abline(0,1, lty=2)
lines(grid, y, lty=2, col=2, type="l")
text(0.9, 0.15, "FDR = 10%", col=2)
abline(0,0.1, lty=1, col=2)

## Conflict detection for split in space:
pLBW.space <- inla.cut(rLBW, split.by="ID.area")

## plot results for split in space
## red line shows false discovery rate (FDR) limit for FDR=10%:
dev.new()
p <- pLBW.space
n.p <- length(p)
grid <- c(0:n.p)/n.p
y <- 0+0.1*grid
sel <- (sort(p) < y[-1])
par(pty="s")
plot(c(1:length(p))/n.p, sort(p), xlim=c(0,1), ylim=c(0,1), pch=19, cex=0.75,
     xlab="Rank of P-values divided by total number (159) of P-values", 
     ylab="Conflict P-values (ordered)", col=c(rep(2, sum(sel)), rep(1, n.p-sum(sel))))
abline(0,1, lty=2)
lines(grid, y, lty=2, col=2, type="l")
text(0.9, 0.15, "FDR = 10%", col=2)
abline(0,0.1, lty=1, col=2)

## Plot map of -log (base 10) of conflict p-values (spatial split).
## Two counties are significantly divergent at FDR=10%
## (Hall and Catoosa, yellow/orange in the map):
dev.new()
g <- georgia
x <- -log(p[order])/log(10)
attr(g, "data")=data.frame(x=x)
spplot(g)

## P-values for the two divergent counties:
print(sort(p)[1:2])
