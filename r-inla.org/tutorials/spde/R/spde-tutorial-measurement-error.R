## ----settings,include=FALSE,message=FALSE, warning=FALSE-----------------
library(knitr)
opts_chunk$set(
fig.path='figs/me',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)
loc.call <- inla.getOption('inla.call')
## inla.setOption(inla.call='remote')
## inla.setOption(num.threads=4)
source('R/spde-tutorial-functions.R')


## ----simloc--------------------------------------------------------------
n.y <- 123;         n.c <- 234
set.seed(1)
loc.c <- cbind(runif(n.c), runif(n.c))
loc.y <- cbind(runif(n.y), runif(n.y))


## ----rfparms-------------------------------------------------------------
kappa.m <- 7;           sigma2.m <- 3
kappa.x <- 10;          sigma2.x <- 2


## ----simula--------------------------------------------------------------
set.seed(2)
mm <- rMatern(n=1, coords=rbind(loc.c, loc.y), 
            kappa=kappa.m, variance=sigma2.m, nu=1)
xx <- rMatern(n=1, coords=loc.y, kappa=kappa.x, 
              variance=sigma2.x, nu=1)
### center it to avoid confounding 
mm <- mm-mean(mm)
xx <- xx-mean(xx)


## ----params--------------------------------------------------------------
alpha.c <- -5;      beta.w <- 0.5
alpha.y <- 3;       beta.c <- 2
sigma2.y <- 0.3


## ----sim-y---------------------------------------------------------------
set.seed(3)
w <- runif(n.c + n.y) 
cc <- alpha.c + beta.w * w + mm 
yy <- alpha.y + beta.c*cc[n.c+1:n.y] + xx + 
    rnorm(n.y, 0, sqrt(sigma2.y))


## ----mesh----------------------------------------------------------------
(rmin <- min(sqrt(8)/c(kappa.m, kappa.x)))
(mesh <- inla.mesh.2d(rbind(loc.c, loc.y), max.edge=rmin/c(5, 2), 
                      cutoff=rmin/10, offset=rmin*c(1/2,3)))$n 


## ----Apred---------------------------------------------------------------
Ac <- inla.spde.make.A(mesh, loc=loc.c)
Ay <- inla.spde.make.A(mesh, loc=loc.y)


## ----dat-----------------------------------------------------------------
stk.y <- inla.stack(data=list(y=cbind(yy, NA, NA)), 
                    A=list(Ay, 1), 
                    effects=list(x=1:mesh$n, 
                        data.frame(
                            a.y=1, o.c=(n.c+1):(n.c+n.y)))) 
stk.c <- inla.stack(data=list(y=cbind(NA, cc[1:n.c], NA)), 
                    A=list(Ac, 1), tag='dat.c', 
                    effects=list(m=1:mesh$n, 
                        data.frame(a.c=1, w=w[1:n.c]))) 
stk.0 <- inla.stack(data=list(y=cbind(NA, NA, rep(0, n.c + n.y))), 
                    A=list(rBind(Ac,Ay), 1), tag='dat.0', 
                    effects=list(m=1:mesh$n, 
                        data.frame(a.c=1, w=w[1:(n.c+n.y)], 
                                   o=1:(n.c+n.y), 
                                   o.weig=rep(-1,n.c+n.y)))) 
stk <- inla.stack(stk.c, stk.y, stk.0) 


## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01


## ----formula-------------------------------------------------------------
form <- y ~  0 + a.c + a.y + w + 
    f(m, model=spde) + f(x, model=spde) + 
    f(o, o.weig, model='iid', 
      hyper=list(theta=list(initial=-20, fixed=TRUE))) + 
    f(o.c, copy='o', fixed=FALSE,  
      hyper=list(theta=list(param=c(0,5)))) 


## ----fit-----------------------------------------------------------------
pcprec <- list(prior='pcprec', param=c(1, 0.01))
res <- inla(form, data=inla.stack.data(stk), family=rep('gaussian',3), 
            control.predictor=list(compute=TRUE, A=inla.stack.A(stk)), 
            control.family=list(list(hyper=list(theta=pcprec)), 
                list(hyper=list(prec=list(initial=20, fixed=TRUE))), 
                list(hyper=list(prec=list(initial=20, fixed=TRUE))))) 


## ----resfix--------------------------------------------------------------
round(cbind(True=c(alpha.c, alpha.y, beta.w), 
            res$summary.fix), 4)


## ----reshyl--------------------------------------------------------------
round(c(True=1/sigma2.y, unlist(res$summary.hy[1,])), 4)


## ----rf------------------------------------------------------------------
round(cbind(True=c(sqrt(8)/kappa.m, sigma2.m, 
                   sqrt(8)/kappa.x, sigma2.x, beta.c), 
            res$summary.hyperpar[-1,]), 3)


## ----likeparsme,eval=F---------------------------------------------------
## par(mfcol=c(2,2), mar=c(3,3,.1,.1), mgp=c(1.5,.5,0), las=1)
## plot(res$marginals.fix[[1]], type='l',
##      xlab=expression(alpha[c]), ylab='')
## abline(v=alpha.c, col=4)
## plot(res$marginals.fix[[2]], type='l',
##      xlab=expression(alpha[y]), ylab='')
## abline(v=alpha.y, col=4)
## plot(res$marginals.fix[[3]], type='l',
##      xlab=expression(beta[w]), ylab='')
## abline(v=beta.w, col=4)
## plot(res$marginals.hy[[6]], type='l',
##      xlab=expression(beta[c]), ylab='')
## abline(v=beta.c, col=4)


## ----vlikeparsme,echo=FALSE,fig.width=7,fig.height=3.5-------------------
par(mfcol=c(2,2), mar=c(3,3,.1,.1), mgp=c(1.5,.5,0), las=1)
plot(res$marginals.fix[[1]], type='l', 
     xlab=expression(alpha[c]), ylab='')
abline(v=alpha.c, col=4)
plot(res$marginals.fix[[2]], type='l', 
     xlab=expression(alpha[y]), ylab='')
abline(v=alpha.y, col=4)
plot(res$marginals.fix[[3]], type='l', 
     xlab=expression(beta[w]), ylab='')
abline(v=beta.w, col=4)
plot(res$marginals.hy[[6]], type='l', 
     xlab=expression(beta[c]), ylab='')
abline(v=beta.c, col=4)


## ----rfparsme,eval=F-----------------------------------------------------
## par(mfcol=c(2,2), mar=c(3,3,.1,.3), mgp=c(1.5,.5,0), las=1)
## for (j in 2:5) {
##     plot(res$marginals.hyperpar[[j]], type='l',
##          xlab=names(res$marginals.hyperpar)[j], ylab='Density')
##     abline(v=c(sqrt(8)/kappa.m, sqrt(sigma2.m),
##                sqrt(8)/kappa.x, sqrt(sigma2.x), beta.c)[j-1], col=4)
## }


## ----vrfparsme,echo=FALSE,fig.width=8,fig.height=4-----------------------
par(mfcol=c(2,2), mar=c(3,3,.1,.3), mgp=c(1.5,.5,0), las=1)
for (j in 2:5) {
    plot(res$marginals.hyperpar[[j]], type='l', 
         xlab=names(res$marginals.hyperpar)[j], ylab='Density')
    abline(v=c(sqrt(8)/kappa.m, sqrt(sigma2.m), 
               sqrt(8)/kappa.x, sqrt(sigma2.x), beta.c)[j-1], col=4)
}


## ----prdmat--------------------------------------------------------------
mesh2locs <- rBind(Ac, Ay)


## ----predicted-----------------------------------------------------------
m.mprd <- drop(mesh2locs%*%res$summary.ran$m$mean)
sd.mprd <- drop(mesh2locs%*%res$summary.ran$m$sd)


## ----prdplotme,eval=F----------------------------------------------------
## plot(m.mprd, mm, asp=1, type='n',
##      xlab='Predicted', ylab='Simulated')
## segments(m.mprd-2*sd.mprd, mm, m.mprd+2*sd.mprd, mm,
##          lty=2, col=gray(.75))
## abline(c(0,1), col=4); points(m.mprd, mm, pch=3, cex=.5)


## ----visprdplotme,echo=FALSE,fig.width=7,fig.height=5--------------------
plot(m.mprd, mm, asp=1, type='n', 
     xlab='Predicted', ylab='Simulated')
segments(m.mprd-2*sd.mprd, mm, m.mprd+2*sd.mprd, mm, 
         lty=2, col=gray(.75))
abline(c(0,1), col=4); points(m.mprd, mm, pch=3, cex=.5)

