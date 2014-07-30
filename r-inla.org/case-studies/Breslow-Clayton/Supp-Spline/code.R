## Section 5.3 B-Spline Nonparametric Regression

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
intKnots <- quantile(unique(bone$age),
                     seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
Omega <- formOmega(a,b,intKnots)
eigOmega <- eigen(Omega)
indsZ <- 1:(numIntKnots+2)
UZ <- eigOmega$vectors[,indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
B <- bs(bone$age,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
ZSpline <- B%*%LZ   
ZBlock <- list(list(groupVec=pdIdent(~ZSpline-1)),list(idnum=pdIdent(~1)))
ZBlock <- unlist(ZBlock,recursive=FALSE)

bone.inla.fit = inla(spnbmd ~ ethnic + age +
        f(idnum,model="iid",param=c(.5, 5e-6),diagonal=0) + 
        f(seqno,model="z",Z=ZSpline,initial=3,param=c(.5, 0.00113)),
        data=bone, family="gaussian", control.predictor=list(compute=TRUE) )
bone.hyperpar = inla.hyperpar (bone.inla.fit)
summary(bone.inla.fit)
inla.emarginal(function(x) 1/x^.5, bone.hyperpar$marginals.hyperpar[[1]])
