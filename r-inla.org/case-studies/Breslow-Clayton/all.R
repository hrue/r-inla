library (INLA)

# file dependency: oxford.RData, salam.RData, spinalbonedata.txt


#################################################################
# Section 5.1 Longitudinal Data - Epilepsy

# visit should be 1 to 4
# Visit is created by glmmAK, it corresponds to Breslow and Clayton's Visit/10, because the codes are -3,-1,1,3. 
require (glmmAK)
data(epilepticBC)
epil = epilepticBC
epil$id2=epil$id
epil$rand=1:nrow(epil)
epil$V4=epil$visit==4
epil$newid=rep(1:(nrow(epil)/4), each=4)

epil.inla.fit.1 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 + f(id,model="iid",param=c(2, 1.140), diagonal=0), data=epil, family="poisson" )
epil.hyperpar.1 = inla.hyperpar(epil.inla.fit.1)
summary(epil.inla.fit.1)
inla.expectation(function(x) 1/x^.5, epil.hyperpar.1$marginals[[1]])

epil.inla.fit.2 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 + f(id,model="iid",param=c(2, 1.240), diagonal=0) +
    f(rand,model="iid",param=c(2, 1.140), diagonal=0), data=epil, family="poisson" )
epil.hyperpar.2 = inla.hyperpar(epil.inla.fit.2)
summary(epil.inla.fit.2)
inla.expectation(function(x) 1/x^.5, epil.hyperpar.2$marginals[[1]])

epil.inla.fit.3 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age + Visit +f(id,  model="2diidwishartpart0", param=c(5, 2.277904, 1.692047, 0), diagonal=0) +
    f(id2, Visit, model="2diidwishartpart1", diagonal = 0), data=epil, family="poisson" )
epil.hyperpar.3 = inla.hyperpar(epil.inla.fit.3)
summary(epil.inla.fit.3)
inla.expectation(function(x) 1/x^.5, epil.hyperpar.3$marginals[[1]])


#################################################################
# Section 5.2 Smoothing of Birth Cohort Effects - Iceland

# dataset downloaded from the winbugs example web page
iceland = list(age = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,
            6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11,
            11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13),
    year = c(6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 10, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7,
                8, 9, 3, 4, 5, 6, 7, 8, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4,
                5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5),
    cases = c(2, 0, 1, 1, 1, 2, 0, 2, 1, 1, 5, 5, 1, 1, 3, 7, 12, 10, 6, 11, 9, 14, 20, 14, 7, 14, 22, 25, 29, 37, 21, 11, 29, 33,
                57, 24, 15, 8, 22, 27, 38, 52, 10, 15, 22, 26, 47, 31, 8, 11, 17, 23, 31, 38, 8, 10, 24, 30, 53, 26, 5, 3, 10, 18,
                22, 30, 1, 7, 11, 26, 32, 17, 5, 8, 17, 32, 31),
    pyr = c(41380, 43650, 49810, 58105, 57105, 76380, 39615, 42205, 48315, 56785, 55965, 33955, 29150, 38460,
                40810, 47490, 55720, 55145, 27950, 37375, 39935, 46895, 54980, 27810, 25055, 27040, 36400, 39355,
            46280, 54350, 24040, 26290, 35480, 38725, 45595, 25710, 22890, 23095, 25410, 34420, 37725, 44740,
            21415, 21870, 24240, 33175, 36345, 21320, 17450, 19765, 20255, 22760, 31695, 34705, 15350, 17720,
            18280, 20850, 29600, 15635, 9965, 12850, 15015, 15725, 18345, 26400, 8175, 11020, 13095, 14050,
            16480, 10885, 7425, 10810, 12260, 14780, 13600)
) 
iceland=data.frame (iceland)
iceland$year2=iceland$year
iceland$year3=iceland$year

x<-cbind(1, 1:11)
extraconstraint<-list(A=matrix((1:11) - mean(1:11), 1, 11), e=matrix(0,1,1))
iceland.inla.fit = inla(cases ~ -1 + as.factor(age) + f(year, model="iid", param=c(0.5, 0.00149) ) + year2 + 
    f(year3, values=1:11, model="rw2", param=c(0.5, 0.001), constr=T, extraconstr=extraconstraint), 
    data=iceland, family="poisson", offset=I(log(pyr)), control.inla= list(int.strategy = "grid",dz=0.05), control.predictor=list(compute=TRUE) ) 
iceland.hyperpar = inla.hyperpar (iceland.inla.fit)
summary(iceland.inla.fit)
inla.expectation(function(x) 1/x^.5, iceland.hyperpar$marginals[[1]])


#################################################################
# Section 5.3 B-Spline Nonparametric Regression

require(splines)
require(nlme)
bone=read.table("spinalbonedata.txt", header=T)
bone=subset(bone, sex=="fem")
bone$ethnic = as.factor(bone$ethnic)
bone$groupVec = 1
bone$seqno=1:nrow(bone)

# Matt Wand code
formOmega <- function(a,b,intKnots)
{
   allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
   K <- length(intKnots) ; L <- 3*(K+8)
   xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
   wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
   Omega     <- t(Bdd*wts)%*%Bdd     
   return(Omega)
}

# Obtain the spline component of the Z matrix:
a <- 8 ; b <- 28; 
numIntKnots <- 15 
intKnots <- quantile(unique(bone$age),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
Omega <- formOmega(a,b,intKnots)
eigOmega <- eigen(Omega)
indsZ <- 1:(numIntKnots+2)
UZ <- eigOmega$vectors[,indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
B <- bs(bone$age,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
ZSpline <- B%*%LZ   
ZBlock <- list(list(groupVec=pdIdent(~ZSpline-1)),list(idnum=pdIdent(~1)))
ZBlock <- unlist(ZBlock,recursive=FALSE)

bone.inla.fit = inla(spnbmd ~ ethnic + age + f(idnum,model="iid",param=c(.5, 5e-6),diagonal=0) + 
    f(seqno,model="z",Z=ZSpline,initial=3,param=c(.5, 0.00113)), data=bone, family="gaussian", control.predictor=list(compute=TRUE) )
bone.hyperpar = inla.hyperpar (bone.inla.fit)
summary(bone.inla.fit)
inla.expectation(function(x) 1/x^.5, bone.hyperpar$marginals[[1]])



#################################################################
# Overdispersion - Deeds

data(Seeds)

seeds.inla.fit.1 = inla(r ~ x1 + x2 + f(plate, model="iid", param=c(.5, .0164)), data=Seeds, family="binomial", Ntrials=n )
seeds.hyperpar.1 = inla.hyperpar(seeds.inla.fit.1)
summary(seeds.inla.fit.1)
inla.expectation(function(x) 1/x^.5, seeds.hyperpar.1$marginals[[1]])

seeds.inla.fit.2 = inla(r ~ x1 + x2 + I(x1*x2) + f(plate, model="iid", param=c(.5, .0164)), data=Seeds, family="binomial", Ntrials=n )
seeds.hyperpar.2 = inla.hyperpar(seeds.inla.fit.2)
summary(seeds.inla.fit.2)
inla.expectation(function(x) 1/x^.5, seeds.hyperpar.2$marginals[[1]])


#################################################################
# Log Odds Ratio - Oxford

load("oxford.Rdata")
oxford$birth.year = oxford$birth.year - 54 # center

oxford.inla.fit = inla( cases ~ as.factor(birth.year) + exposed + I(exposed * birth.year) + I(exposed * (birth.year**2-22)) +
    f(birth.year, model="iid", param=c(.5, .0164), weights=exposed), data=oxford, family="binomial", Ntrials=total ) 
oxford.hyperpar = inla.hyperpar (oxford.inla.fit)
summary(oxford.inla.fit)
inla.expectation(function(x) 1/x^.5, oxford.hyperpar$marginals[[1]])


#################################################################
# Spatial Aggregation - Scotland

data(Scotland)
Scotland$Region2=Scotland$Region

scotland.inla.fit.1 = inla(Counts ~ 1 + f(Region, model="iid", param=c(1, .0014)), data=Scotland, family="poisson", E=E )
scotland.hyperpar.1 = inla.hyperpar (scotland.inla.fit.1)
summary(scotland.inla.fit.1)
inla.expectation(function(x) 1/x^.5, scotland.hyperpar.1$marginals[[1]])

scotland.inla.fit.2 = inla(Counts ~ 1 + I(X/10) + f(Region, model="iid", param=c(1, .0014)), data=Scotland, family="poisson", E=E )
scotland.hyperpar.2 = inla.hyperpar (scotland.inla.fit.2)
summary(scotland.inla.fit.2)
inla.expectation(function(x) 1/x^.5, scotland.hyperpar.2$marginals[[1]])

scotland.inla.fit.3 = inla(Counts ~ 1 + f(Region, model="iid", param=c(1, .0014)) +
    f(Region2, model="besag", graph="scotland.graph", param=c(1, .2/.59)), data=Scotland, family="poisson", E=E )
scotland.hyperpar.3 = inla.hyperpar (scotland.inla.fit.3)
summary(scotland.inla.fit.3)
inla.expectation(function(x) 1/x^.5, scotland.hyperpar.3$marginals[[1]])

scotland.inla.fit.4 = inla(Counts ~ 1 + I(X/10) + f(Region, model="iid", param=c(1, .0014)) +
    f(Region2, model="besag", graph="scotland.graph", param=c(1, .4/.59)), data=Scotland, family="poisson", E=E )
scotland.hyperpar.4 = inla.hyperpar (scotland.inla.fit.4)
summary(scotland.inla.fit.4)
inla.expectation(function(x) 1/x^.5, scotland.hyperpar.4$marginals[[1]])


#################################################################
# Crossed Random Effects - Salamander

load("salam.RData")
# organize data into a form suitable for logistic regression
dat0=data.frame("y"=c(salam$y), "fW"=as.integer(salam$x[,"W/R"]==1 | salam$x[,"W/W"]==1), 
    "mW"=as.integer(salam$x[,"R/W"]==1 | salam$x[,"W/W"]==1), 
    "WW"=as.integer(salam$x[,"W/W"]==1 ) )
# add salamander id
id = t( apply(salam$z, 1, function(x) {
        tmp = which (x==1)
        tmp[2] = tmp[2] - 20
        tmp
    }) ) 
# ids are suitable for model A and C, but not B
id.modA = rbind(id, id+40, id+20)
colnames (id.modA) = c("f.modA","m.modA")
dat0=cbind (dat0, id.modA, group=1)
dat0$experiment=as.factor(rep(1:3, each=120))
dat0$group=as.factor(dat0$group)

salamander = dat0
salamander.e1 = subset (dat0, dat0$experiment==1)
salamander.e2 = subset (dat0, dat0$experiment==2)
salamander.e3 = subset (dat0, dat0$experiment==3)

# salamander.e1
salamander.e1.inla.fit = inla(y~fW+mW+WW + f(f.modA, model="iid", param=c(1,.622)) + f(m.modA, model="iid", param=c(1,.622)), 
    family="binomial", data=salamander.e1, Ntrials=rep(1,nrow(salamander.e1)))
salamander.e1.hyperpar = inla.hyperpar (salamander.e1.inla.fit)
summary(salamander.e1.inla.fit)
summary(salamander.e1.hyperpar)

inla.expectation(function(x) 1/x^.5, salamander.e1.hyperpar$marginals[[1]])
inla.expectation(function(x) 1/x^.5, salamander.e1.hyperpar$marginals[[2]])

# same for salamander.e2 and salamander.e3
