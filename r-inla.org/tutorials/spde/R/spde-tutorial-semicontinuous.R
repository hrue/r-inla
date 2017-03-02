### R code from vignette source 'spde-tutorial-semicontinuous.Rnw'

###################################################
### code chunk number 1: sett
###################################################
options(width=75, prompt = " ", continue = "   ") 
library(splancs) 
library(INLA) 
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)
library(gridExtra) 
library(lattice) 


###################################################
### code chunk number 2: readdata
###################################################
data(PRprec) 


###################################################
### code chunk number 3: prechead
###################################################
PRprec[1:3, 1:10] 


###################################################
### code chunk number 4: summ
###################################################
sapply(PRprec[,4:10], summary) 
colSums(PRprec[,4:10]>0, na.rm=TRUE) 


###################################################
### code chunk number 5: vardef
###################################################
jday <- 6 
z <- (PRprec[,jday]>0) + 0 
y <- ifelse(PRprec[,jday]>0, PRprec[,jday], NA) 


###################################################
### code chunk number 6: rain (eval = FALSE)
###################################################
## par(mfrow=c(1, 2), mar=c(3,3,.5,0), mgp=c(1.5,0.5,0)) 
## hist(y, xlab='mm', col=gray(.9), main='') 
## par(mar=c(0,0,0,0)) 
## plot(PRprec[,1:2], cex=0.3 + PRprec[,jday]/20, asp=1, 
##      col=c('red', 'blue')[z+1], axes=FALSE) 
## q.y <- c(0, quantile(y, 0:7/7, na.rm=T)) 
## legend('topright', format(q.y, dig=2), pch=1, bty='n', 
##        pt.cex=0.3+q.y/20, col=c('red',rep('blue',length(q.y)-1))) 


###################################################
### code chunk number 7: vrain
###################################################
par(mfrow=c(1, 2), mar=c(3,3,.5,0), mgp=c(1.5,0.5,0)) 
hist(y, xlab='mm', col=gray(.9), main='') 
par(mar=c(0,0,0,0)) 
plot(PRprec[,1:2], cex=0.3 + PRprec[,jday]/20, asp=1, 
     col=c('red', 'blue')[z+1], axes=FALSE) 
q.y <- c(0, quantile(y, 0:7/7, na.rm=T)) 
legend('topright', format(q.y, dig=2), pch=1, bty='n', 
       pt.cex=0.3+q.y/20, col=c('red',rep('blue',length(q.y)-1))) 


###################################################
### code chunk number 8: spde-tutorial-semicontinuous.Rnw:181-184
###################################################
mesh <- inla.mesh.2d(cbind(PRprec[,1], PRprec[,2]), 
                     max.edge=c(0.5, 1), cutoff=0.05) 
mesh$n


###################################################
### code chunk number 9: spdedef
###################################################
spde <- inla.spde2.matern(mesh, alpha=2) 


###################################################
### code chunk number 10: A
###################################################
A <- inla.spde.make.A(mesh, loc=as.matrix(PRprec[,1:2])) 


###################################################
### code chunk number 11: modelssep
###################################################
stk.z <- inla.stack(tag='est.z', 
                    data=list(z=z, ### occurrence for separate model
                        y=cbind(z, NA)), ### z at first column of y
                    A=list(A, 1), 
                    effects=list( 
                      list(i.z=1:spde$n.spde), 
                      list(z.b0=rep(1,length(z))))) 


###################################################
### code chunk number 12: stack
###################################################
stk.y <- inla.stack(tag='est.y', 
                    data=list(r=y, ### rainfall for separate model
                        y=cbind(NA, y)), ### rainfall at second column 
                    A=list(A, 1), 
                    effects=list( 
                      list(i.y=1:spde$n.spde), 
                      list(y.b0=rep(1,length(y))))) 


###################################################
### code chunk number 13: ressep
###################################################
res.z <- inla(z ~ 0 + z.b0 + f(i.z, model=spde), family='binomial', 
              data=inla.stack.data(stk.z), control.compute=list(dic=TRUE), 
              control.predictor=list(A=inla.stack.A(stk.z), compute=TRUE)) 
res.y <- inla(r ~ 0 +  y.b0 + f(i.y, model=spde), family='gamma', 
              data=inla.stack.data(stk.y), control.compute=list(dic=TRUE), 
              control.predictor=list(A=inla.stack.A(stk.y), compute=TRUE)) 


###################################################
### code chunk number 14: stk.full
###################################################
stk.zy <- inla.stack(stk.z, stk.y) 


###################################################
### code chunk number 15: model
###################################################
res.zy <- inla(y ~ 0 + z.b0 + y.b0 + 
               f(i.z, model=spde) + f(i.y, copy='i.z', fixed=FALSE), 
               family=c('binomial', 'gamma'), 
               data=inla.stack.data(stk.zy), control.compute=list(dic=TRUE), 
               control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE)) 


###################################################
### code chunk number 16: modelyz
###################################################
res.yz <- inla(y ~ 0 + z.b0 + y.b0 + 
               f(i.y, model=spde) + f(i.z, copy='i.y', fixed=FALSE), 
               family=c('binomial', 'gamma'), 
               data=inla.stack.data(stk.zy), control.compute=list(dic=TRUE), 
               control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE)) 


###################################################
### code chunk number 17: beta-spde
###################################################
round(rbind(beta.zy=res.zy$summary.hy[4, ], beta.yz=res.yz$summary.hy[4, ]), 4)


###################################################
### code chunk number 18: dics
###################################################
iz <- inla.stack.index(stk.zy, tag='est.z')$data 
dic.zy <- c(dic.z = sum(res.zy$dic$local.dic[iz], na.rm=TRUE), 
            dic.y = sum(res.zy$dic$local.dic[-iz], na.rm=TRUE))
dic.yz <- c(dic.z = sum(res.yz$dic$local.dic[iz], na.rm=TRUE), 
            dic.y = sum(res.yz$dic$local.dic[-iz], na.rm=TRUE))


###################################################
### code chunk number 19: dics-sj
###################################################
rbind(sep=c(res.z$dic$dic, res.y$dic$dic), joint.zy=dic.zy, joint.yz=dic.yz)     


###################################################
### code chunk number 20: summaries-alphaz
###################################################
round(rbind(sep=res.z$summary.fix, joint.zy=res.zy$summary.fix[1,], 
            joint.yz=res.yz$summary.fix[1,]), 4)


###################################################
### code chunk number 21: p0alphaz
###################################################
sapply(list(sep=res.z$marginals.fix[[1]], joint.zy=res.zy$marginals.fix[[1]], 
            joint.yz=res.yz$marginals.fix[[1]]), function(m) inla.pmarginal(0, m))


###################################################
### code chunk number 22: summaries-alphaz
###################################################
round(rbind(sep=res.y$summary.fix, joint.zy=res.zy$summary.fix[2,], 
            joint.yz=res.yz$summary.fix[2,]), 4)


###################################################
### code chunk number 23: invlink
###################################################
c(binomial(link='logit')$linkinv(res.zy$summary.fix[1,1]), 
  exp(res.zy$summary.fix[2,1])) 


###################################################
### code chunk number 24: obs
###################################################
c(occ=mean(z, na.rm=TRUE), rain=mean(y, na.rm=TRUE)) 


###################################################
### code chunk number 25: sphi
###################################################
res.yz$summary.hy[1, ] 


###################################################
### code chunk number 26: post (eval = FALSE)
###################################################
## par(mfrow=c(2,2), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0), las=1) 
## plot(res.yz$marginals.fix[[1]], type='l', ylab='Density', 
##      xlab=expression(alpha[z])) 
## plot(res.yz$marginals.fix[[2]], type='l', ylab='Density', 
##      xlab=expression(alpha[y])) 
## plot(res.yz$marginals.hy[[1]], type='l', ylab='Density', 
##      xlab=expression(phi)) 
## plot(res.yz$marginals.hy[[4]], type='l', ylab='Density', 
##      xlab=expression(beta)) 


###################################################
### code chunk number 27: vpost
###################################################
par(mfrow=c(2,2), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0), las=1) 
plot(res.yz$marginals.fix[[1]], type='l', ylab='Density', 
     xlab=expression(alpha[z])) 
plot(res.yz$marginals.fix[[2]], type='l', ylab='Density', 
     xlab=expression(alpha[y])) 
plot(res.yz$marginals.hy[[1]], type='l', ylab='Density', 
     xlab=expression(phi)) 
plot(res.yz$marginals.hy[[4]], type='l', ylab='Density', 
     xlab=expression(beta)) 


###################################################
### code chunk number 28: xq (eval = FALSE)
###################################################
## ordx <- order(res.yz$summary.random$i.z$mean) 
## par(mar=c(3,3,0.5,0.5), mgp=c(1.5,0.5,0), las=1) 
## plot(res.yz$summary.random$i.z$mean[ordx], type='l', ylab='x', 
##      ylim=range(res.yz$summary.random$i.z[, 4:6])) 
## for (i in c(4,6)) lines(res.yz$summary.random$i.z[ordx,i], lty=2) 
## abline(h=0, lty=3) 


###################################################
### code chunk number 29: vxq
###################################################
ordx <- order(res.yz$summary.random$i.z$mean) 
par(mar=c(3,3,0.5,0.5), mgp=c(1.5,0.5,0), las=1) 
plot(res.yz$summary.random$i.z$mean[ordx], type='l', ylab='x', 
     ylim=range(res.yz$summary.random$i.z[, 4:6])) 
for (i in c(4,6)) lines(res.yz$summary.random$i.z[ordx,i], lty=2) 
abline(h=0, lty=3) 


###################################################
### code chunk number 30: refit
###################################################
res.yz0 <- inla(y ~ 0 + z.b0 + y.b0, 
            family=c('binomial', 'gamma'), 
            data=inla.stack.data(stk.zy), control.compute=list(dic=TRUE), 
            control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE)) 


###################################################
### code chunk number 31: dics
###################################################
rbind(dic.0=c(dic.yz0=sum(res.yz0$dic$local.dic[iz], na.rm=TRUE), 
          dic.yz0=sum(res.yz0$dic$local.dic[-iz], na.rm=TRUE)), dic.s=dic.yz)


###################################################
### code chunk number 32: resfield
###################################################
res.yz.f <- inla.spde2.result(res.yz, 'i.y', spde, do.transf=TRUE) 


###################################################
### code chunk number 33: mvar
###################################################
inla.emarginal(function(x) x, res.yz.f$marginals.variance.nominal[[1]]) 


###################################################
### code chunk number 34: range
###################################################
inla.emarginal(function(x) x, res.yz.f$marginals.range.nominal[[1]]) 


###################################################
### code chunk number 35: xpars (eval = FALSE)
###################################################
## par(mfrow=c(1,3), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0) 
## plot.default(inla.tmarginal(function(x) exp(x), res.yz$marginals.hy[[3]]), 
##              type='l', xlab=expression(kappa), ylab='Density') 
## plot.default(res.yz.f$marginals.variance.nominal[[1]], type='l', 
##              xlab=expression(sigma[x]^2), ylab='Density') 
## plot.default(res.yz.f$marginals.range.nominal[[1]], type='l', 
##              xlab='Practical range', ylab='Density') 


###################################################
### code chunk number 36: vxpars
###################################################
par(mfrow=c(1,3), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0) 
plot.default(inla.tmarginal(function(x) exp(x), res.yz$marginals.hy[[3]]), 
             type='l', xlab=expression(kappa), ylab='Density') 
plot.default(res.yz.f$marginals.variance.nominal[[1]], type='l', 
             xlab=expression(sigma[x]^2), ylab='Density') 
plot.default(res.yz.f$marginals.range.nominal[[1]], type='l', 
             xlab='Practical range', ylab='Density') 


###################################################
### code chunk number 37: projgrid
###################################################
data(PRborder) 
(nxy <- round(c(diff(range(PRborder[,1])), diff(range(PRborder[,2])))/.02)) 
projgrid <- inla.mesh.projector(mesh, xlim=range(PRborder[,1]), 
                                ylim=range(PRborder[,2]), dims=nxy) 


###################################################
### code chunk number 38: projpred
###################################################
xmean <- inla.mesh.project(projgrid, res.yz$summary.random$i.z$mean) 
xsd <- inla.mesh.project(projgrid, res.yz$summary.random$i.z$sd) 


###################################################
### code chunk number 39: sp
###################################################
table(xy.in <- inout(projgrid$lattice$loc, PRborder)) 
xmean[!xy.in] <- xsd[!xy.in] <- NA 


###################################################
### code chunk number 40: xrain3c (eval = FALSE)
###################################################
## grid.arrange(levelplot(xmean, col.regions=topo.colors(99), 
##                        xlab='', ylab='', scales=list(draw=FALSE)), 
##              levelplot(xsd, col.regions=topo.colors(99), 
##                        xlab='', ylab='', scales=list(draw=FALSE)), nrow=1) 


###################################################
### code chunk number 41: xrain3
###################################################
grid.arrange(levelplot(xmean, col.regions=topo.colors(99), 
                       xlab='', ylab='', scales=list(draw=FALSE)), 
             levelplot(xsd, col.regions=topo.colors(99), 
                       xlab='', ylab='', scales=list(draw=FALSE)), nrow=1) 


