### R code from vignette source 'spde-tutorial-jcovar.Rnw'

###################################################
### code chunk number 1: sett
###################################################
options(width=75, prompt = " ", continue = "   ")
library(INLA)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)
library(geoR)


###################################################
### code chunk number 2: params
###################################################
n.y <- 123;         n.c <- 234
alpha.y <- -5;      beta <- -3
sigma2.y <- 0.5;    sigma2.c <- 0.2


###################################################
### code chunk number 3: simloc
###################################################
set.seed(1)
loc.c <- cbind(runif(n.c), runif(n.c))
loc.y <- cbind(runif(n.y), runif(n.y))


###################################################
### code chunk number 4: rfparms
###################################################
kappa.m <- 3;          sigma2.m <- 5
kappa.x <- 7;          sigma2.x <- 3


###################################################
### code chunk number 5: simula
###################################################
library(geoR)
set.seed(2)
mm <- grf(grid=rbind(loc.c, loc.y), messages=FALSE, 
          cov.pars=c(sigma2.m, 1/kappa.m), kappa=1)$data
xx <- grf(grid=loc.y, messages=FALSE, 
          cov.pars=c(sigma2.x, 1/kappa.x), kappa=1)$data


###################################################
### code chunk number 6: sim-y
###################################################
set.seed(3)
cc <- mm + rnorm(n.c+n.y, 0, sqrt(sigma2.c))
yy <- alpha.y + beta*mm[n.c+1:n.y] + xx + 
  rnorm(n.y, 0, sqrt(sigma2.y))


###################################################
### code chunk number 7: mesh
###################################################
pl01 <- matrix(c(0,1,1,0,0, 0,0,1,1,0), ncol=2)
(mesh <- inla.mesh.create.helper(, pl01, cutoff=.05, 
                                 max.edge=c(.1,.2)))$n


###################################################
### code chunk number 8: Apred
###################################################
Ac <- inla.spde.make.A(mesh, loc=loc.c)
Ay <- inla.spde.make.A(mesh, loc=loc.y)


###################################################
### code chunk number 9: dat
###################################################
stk.c <- inla.stack(data=list(y=cbind(cc[1:n.c], NA)), 
                    A=list(Ac), tag='dat.c', 
                    effects=list(m=1:mesh$n))
stk.y <- inla.stack(data=list(y=cbind(NA, yy)), 
                    A=list(Ay, 1), 
                    effects=list(list(c.y=1:mesh$n, x=1:mesh$n), 
                      list(a.y=rep(1,length(yy)))))
stk <- inla.stack(stk.c, stk.y)


###################################################
### code chunk number 10: formula
###################################################
form <- y ~  0 + f(m, model=spde) + 
  a.y + f(x, model=spde) + f(c.y, copy='m', fixed=FALSE, 
               hyper=list(theta=list(param=c(-3,25), initial=-3))) 


###################################################
### code chunk number 11: fit
###################################################
spde <- inla.spde2.matern(mesh)
res <- inla(form, data=inla.stack.data(stk), family=rep('gaussian',2), 
            control.predictor=list(compute=TRUE, A=inla.stack.A(stk))) 


###################################################
### code chunk number 12: resfix
###################################################
round(cbind(True=alpha.y, res$summary.fix), 4)


###################################################
### code chunk number 13: reshyl
###################################################
round(cbind(True=1/c(Prec.c=sigma2.c, Prec.y=sigma2.y), 
            res$summary.hy[1:2,]), 4)


###################################################
### code chunk number 14: reshyb
###################################################
round(cbind(True=beta, res$summary.hy[7,]), 4)


###################################################
### code chunk number 15: rf
###################################################
m.rf <- inla.spde2.result(res, 'm', spde)
x.rf <- inla.spde2.result(res, 'x', spde)
round(cbind(True=c(s2m=sigma2.m, s2x=sigma2.x), 
            mean=c(inla.emarginal(function(x) x, 
              m.rf$marginals.variance.n[[1]]), 
              inla.emarginal(function(x) x, 
                             x.rf$marginals.variance.n[[1]])), 
            rbind(inla.hpdmarginal(.95, m.rf$marginals.variance.n[[1]]), 
                  inla.hpdmarginal(.95, x.rf$marginals.variance.n[[1]]))), 4)


###################################################
### code chunk number 16: kappa
###################################################
post.k <- lapply(res$marginals.hy[c(4,6)], function(y) 
                 inla.tmarginal(function(x) exp(x), y))
round(cbind(True=c(kappa.m=kappa.m, kappa.x=kappa.x), 
            t(sapply(post.k, function(x) 
                     c(mean=inla.emarginal(function(y) y, x), 
                       CI=inla.hpdmarginal(.95, x))))), 4)


###################################################
### code chunk number 17: likeparsjcov (eval = FALSE)
###################################################
## par(mfcol=c(2,2), mar=c(3,3,.1,.1), mgp=c(1.5,.5,0), las=1)
## plot(res$marginals.fix[[1]], type='l', 
##      xlab=expression(alpha[y]), ylab='')
## abline(v=alpha.y, col=4)
## plot(res$marginals.hy[[7]], type='l', 
##      xlab=expression(beta), ylab='')
## abline(v=beta, col=4)
## plot.default(inla.tmarginal(function(x) 1/x, res$marginals.hy[[1]]), 
##              type='l', xlab=expression(sigma[c]^2), ylab='')
## abline(v=sigma2.c, col=4)
## plot.default(inla.tmarginal(function(x) 1/x, res$marginals.hy[[2]]), 
##              type='l', xlab=expression(sigma[y]^2), ylab='')
## abline(v=sigma2.y, col=4)


###################################################
### code chunk number 18: likeparsjcov
###################################################
par(mfcol=c(2,2), mar=c(3,3,.1,.1), mgp=c(1.5,.5,0), las=1)
plot(res$marginals.fix[[1]], type='l', 
     xlab=expression(alpha[y]), ylab='')
abline(v=alpha.y, col=4)
plot(res$marginals.hy[[7]], type='l', 
     xlab=expression(beta), ylab='')
abline(v=beta, col=4)
plot.default(inla.tmarginal(function(x) 1/x, res$marginals.hy[[1]]), 
             type='l', xlab=expression(sigma[c]^2), ylab='')
abline(v=sigma2.c, col=4)
plot.default(inla.tmarginal(function(x) 1/x, res$marginals.hy[[2]]), 
             type='l', xlab=expression(sigma[y]^2), ylab='')
abline(v=sigma2.y, col=4)


###################################################
### code chunk number 19: rfparsjcov (eval = FALSE)
###################################################
## par(mfcol=c(2,3), mar=c(3,3,.1,.3), mgp=c(1.5,.5,0), las=1)
## plot.default(post.k[[1]], type='l', xlab=expression(kappa[m]), ylab='')
## abline(v=kappa.m, col=4)
## plot.default(post.k[[2]], type='l', xlab=expression(kappa[x]), ylab='')
## abline(v=kappa.x, col=4)
## plot.default(m.rf$marginals.variance.n[[1]], type='l', 
##              xlab=expression(sigma[m]^2), ylab='')
## abline(v=sigma2.m, col=4)
## plot.default(x.rf$marginals.variance.n[[1]], type='l', 
##              xlab=expression(sigma[x]^2), ylab='')
## abline(v=sigma2.x, col=4)
## plot.default(m.rf$marginals.range.n[[1]], type='l', 
##              xlab='practical range of m', ylab='')
## abline(v=sqrt(8)/kappa.m, col=4)
## plot.default(x.rf$marginals.range.n[[1]], type='l', 
##              xlab='practical range of x', ylab='')
## abline(v=sqrt(8)/kappa.x, col=4)


###################################################
### code chunk number 20: vrfparsjcov
###################################################
par(mfcol=c(2,3), mar=c(3,3,.1,.3), mgp=c(1.5,.5,0), las=1)
plot.default(post.k[[1]], type='l', xlab=expression(kappa[m]), ylab='')
abline(v=kappa.m, col=4)
plot.default(post.k[[2]], type='l', xlab=expression(kappa[x]), ylab='')
abline(v=kappa.x, col=4)
plot.default(m.rf$marginals.variance.n[[1]], type='l', 
             xlab=expression(sigma[m]^2), ylab='')
abline(v=sigma2.m, col=4)
plot.default(x.rf$marginals.variance.n[[1]], type='l', 
             xlab=expression(sigma[x]^2), ylab='')
abline(v=sigma2.x, col=4)
plot.default(m.rf$marginals.range.n[[1]], type='l', 
             xlab='practical range of m', ylab='')
abline(v=sqrt(8)/kappa.m, col=4)
plot.default(x.rf$marginals.range.n[[1]], type='l', 
             xlab='practical range of x', ylab='')
abline(v=sqrt(8)/kappa.x, col=4)


###################################################
### code chunk number 21: prdmat
###################################################
mesh2locs <- rBind(Ac, Ay)


###################################################
### code chunk number 22: predicted
###################################################
m.mprd <- drop(mesh2locs%*%res$summary.ran$m$mean)
sd.mprd <- drop(mesh2locs%*%res$summary.ran$m$sd)


###################################################
### code chunk number 23: mprdplot (eval = FALSE)
###################################################
## plot(m.mprd, mm, asp=1, type='n', 
##      xlab='Predicted', ylab='Simulated')
## segments(m.mprd-2*sd.mprd, mm, m.mprd+2*sd.mprd, mm, 
##          lty=2, col=gray(.75))
## abline(c(0,1), col=4); points(m.mprd, mm, pch=3, cex=.5)


###################################################
### code chunk number 24: vismprdplot
###################################################
plot(m.mprd, mm, asp=1, type='n', 
     xlab='Predicted', ylab='Simulated')
segments(m.mprd-2*sd.mprd, mm, m.mprd+2*sd.mprd, mm, 
         lty=2, col=gray(.75))
abline(c(0,1), col=4); points(m.mprd, mm, pch=3, cex=.5)


