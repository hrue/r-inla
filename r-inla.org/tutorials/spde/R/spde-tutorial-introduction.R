## ----sett,echo=FALSE,results='hide',message=FALSE,warning=FALSE----------
library(knitr)
opts_chunk$set(
fig.path='figs/intro',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)
source('R/spde-tutorial-functions.R')

## ----maternsamplefunc,results='hide'-------------------------------------
cMatern <- function(h, nu, kappa) ### Matern correlation
    besselK(h * kappa, nu) *
        (h*kappa)^nu / (gamma(nu) * 2^(nu-1))
### function to sample from zero mean multivariate normal
rmvnorm0 <- function(n, cov, L=NULL) { 
    if (is.null(L)) L <- chol(cov)
    return(crossprod(L, matrix(rnorm(n*ncol(L)), ncol(L))))
}

## ----loc1----------------------------------------------------------------
### define locations and distance matrix
loc <- 1:249/25 
mdist <- as.matrix(dist(loc))

## ----param---------------------------------------------------------------
### define parameters
nu <- c(0.5,1,2,5)
pract.range = c(1,4) 
kappa <- c(sqrt(8*nu)/pract.range[1], 
           sqrt(8*nu)/pract.range[2]) 
### covariance parameter scenarios
params <- cbind(nu=rep(nu, length(pract.range)),
                kappa=kappa,
                r=rep(pract.range, each=length(nu)))

## ----error---------------------------------------------------------------
### sample error
set.seed(123)
z <- matrix(rnorm(nrow(mdist)*5), ncol=5)

## ----samples-------------------------------------------------------------
### compute the correlated samples
yy <- lapply(1:nrow(params), function(j) { ## scenarios
    v <- cMatern(mdist, params[j,1], params[j,2])
    diag(v) <- 1 + 1e-10
    return(list(params=params[j,], ### parameter scenario
                y=crossprod(chol(v), z))) ### compute sample
})

## ----maternsamples,eval=FALSE,results='hide'-----------------------------
## ### visualize
## (ry <- range(unlist(lapply(yy, tail, 1))))
## par(mfcol=c(4,2), mar=c(2,2,1,.1), mgp=c(1.5,0.7,0), las=1)
## for (i in 1:length(yy)) { ### each scenario
##     plot(loc, yy[[i]]$y[,1], ylim=ry,
##          xlab='', ylab='', type='n',
##          main=as.expression(bquote(paste(
##              nu==.(yy[[i]]$params[1]), ', ',
##              kappa==.(round(yy[[i]]$params[2],2)), ', ',
##              r==.(yy[[i]]$params[3])))))
##     for (k in 1:5)
##         lines(loc, yy[[i]]$y[,k], col=k) ### each sample
## }

## ----maternsamplesfig,echo=FALSE,results='hide',fig.width=5,fig.height=7,out.width='0.97\\textwidth'----
### visualize
(ry <- range(unlist(lapply(yy, tail, 1))))
par(mfcol=c(4,2), mar=c(2,2,1,.1), mgp=c(1.5,0.7,0), las=1)
for (i in 1:length(yy)) { ### each scenario
    plot(loc, yy[[i]]$y[,1], ylim=ry, 
         xlab='', ylab='', type='n',
         main=as.expression(bquote(paste(
             nu==.(yy[[i]]$params[1]), ', ',
             kappa==.(round(yy[[i]]$params[2],2)), ', ',
             r==.(yy[[i]]$params[3])))))
    for (k in 1:5) 
        lines(loc, yy[[i]]$y[,k], col=k) ### each sample 
}

## ----rpts----------------------------------------------------------------
n <- 200;  set.seed(123) 
pts <- cbind(s1=sample(1:n/n-0.5/n)^2, s2=sample(1:n/n-0.5/n)^2)

## ----distpts-------------------------------------------------------------
dmat <- dist(pts)

## ----params--------------------------------------------------------------
beta0 <- 10; sigma2e <- 0.3; sigma2u <- 5; kappa <- 7; nu <- 1

## ----covMatm-------------------------------------------------------------
mcor <- as.matrix(2^(1-nu)*(kappa*dmat)^nu * 
                  besselK(dmat*kappa,nu)/gamma(nu)) 
diag(mcor) <- 1;   mcov <- sigma2e*diag(n) + sigma2u*mcor 

## ----chol1mvnorm---------------------------------------------------------
L <- chol(mcov);   set.seed(234) 
y1 <- beta0 + drop(crossprod(L, rnorm(n))) 

## ----plot1c,eval=FALSE---------------------------------------------------
## par(mar=c(3,3,1,1), mgp=c(1.7, 0.7, 0), las=1)
## plot(pts, asp=1, xlim=c(0,1.2), cex=y1/10)
## q <- quantile(y1, 0:5/5)
## legend('topright', format(q, dig=2), pch=1, pt.cex=q/10)

## ----plot1,echo=FALSE,fig.width=5.5,fig.height=4.7-----------------------
par(mar=c(3,3,1,1), mgp=c(1.7, 0.7, 0), las=1)
plot(pts, asp=1, xlim=c(0,1.2), cex=y1/10)
q <- quantile(y1, 0:5/5)
legend('topright', format(q, dig=2), pch=1, pt.cex=q/10)

## ----datatoy-------------------------------------------------------------
data(SPDEtoy)

## ----rw1rw2--------------------------------------------------------------
(q1 <- INLA:::inla.rw1(n=5))
crossprod(q1) ### same inner pattern as for RW2
INLA:::inla.rw2(n=5)

## ----mesh0, echo=TRUE, results='hide'------------------------------------
s <- 3 ### this factor will only changes C, not G
pts <- rbind(c(1,1), c(2,1), 
             c(2.6, 1), c(0.7,1.7), 4:5/3, c(2,1.7))*s
n <- nrow(pts)
mesh0 <- inla.mesh.2d(pts[1:2,], max.edge=3*s, 
                      offset=1*s, n=6, cutoff=s*1/2)
mesh <- inla.mesh.2d(rbind(c(3.3,1)*s, c(2.4,2)*s, 
                           mesh0$loc[-c(3:4),1:2]), 
                     max.edge=3*s, offset=1e-5, cutoff=s*1/2, n=100)
(m <- mesh$n)
dmesh <- inla.mesh.dual(mesh)
fem <- inla.mesh.fem(mesh, order=1)
A <- inla.spde.make.A(mesh, pts)

## ----mesh0plot,eval=FALSE------------------------------------------------
## par(mfrow=c(1,3), mar=c(2,2,1,1))
## plot(mesh, asp=1, lwd=2, edge.color=1)
## box(); axis(1); axis(2)
## points(mesh$loc, cex=3)
## text(mesh$loc[,1]-rep(c(0,0.1),c(m-1,1))*s,
##      mesh$loc[,2]+.2*s, 1:m, cex=2)
## 
## plot(dmesh, asp=1, lwd=2, main='Dual mesh overlayed')
## plot(mesh, add=TRUE)
## box(); axis(1); axis(2)
## points(mesh$loc, cex=3)
## 
## plot(mesh, asp=1, lwd=2, edge.color=1, main='')
## title(main='Mesh and points')
## box(); axis(1); axis(2)
## points(mesh$loc, cex=3)
## points(pts, pch=8, cex=2, lwd=2)
## text(pts[,1], pts[,2]+0.2*s, 1:n, cex=2)
## 
## A <- as(Matrix(round(as.matrix(A), 10)), ### force zero as zero
##         'dgTMatrix')
## cc <- as(fem$c0, 'dgTMatrix')
## gg <- as(fem$g1, 'dgTMatrix')
## library(gridExtra)
## library(latticeExtra)
## grid.arrange(
##     plot(cc, colorkey=FALSE, xlab='', ylab='', sub='') +
##     layer(panel.text(cc@j+1L, cc@i+1L, paste0(round(cc@x)),
##                      col=gray(cc@x>30))),
##     plot(gg, colorkey=FALSE, xlab='', ylab='', sub='') +
##     layer(panel.text(gg@j+1L, gg@i+1L, round(gg@x,2), col=gray(abs(gg@x)>1.5))),
##     plot(A, colorkey=FALSE, xlab='', ylab='', sub='') +
##     layer(panel.text(A@j+1L, A@i+1L, round(A@x, 2),
##                      col=gray(A@x>0.5))), ncol=3)

## ----mesh0fig,echo=FALSE,results='hide',fig.width=12,fig.height=4, out.width='0.99\\linewidth'----
par(mfrow=c(1,3), mar=c(2,2,1,1))
plot(mesh, asp=1, lwd=2, edge.color=1)
box(); axis(1); axis(2)
points(mesh$loc, cex=3)
text(mesh$loc[,1]-rep(c(0,0.1),c(m-1,1))*s, 
     mesh$loc[,2]+.2*s, 1:m, cex=2)

plot(dmesh, asp=1, lwd=2, main='Dual mesh overlayed')
plot(mesh, add=TRUE)
box(); axis(1); axis(2)
points(mesh$loc, cex=3)

plot(mesh, asp=1, lwd=2, edge.color=1, main='')
title(main='Mesh and points')
box(); axis(1); axis(2)
points(mesh$loc, cex=3)
points(pts, pch=8, cex=2, lwd=2)
text(pts[,1], pts[,2]+0.2*s, 1:n, cex=2)

A <- as(Matrix(round(as.matrix(A), 10)), ### force zero as zero
        'dgTMatrix')
cc <- as(fem$c0, 'dgTMatrix')
gg <- as(fem$g1, 'dgTMatrix')
library(gridExtra)
library(latticeExtra)
grid.arrange(
    plot(cc, colorkey=FALSE, xlab='', ylab='', sub='') + 
    layer(panel.text(cc@j+1L, cc@i+1L, paste0(round(cc@x)), 
                     col=gray(cc@x>30))), 
    plot(gg, colorkey=FALSE, xlab='', ylab='', sub='') + 
    layer(panel.text(gg@j+1L, gg@i+1L, round(gg@x,2), col=gray(abs(gg@x)>1.5))),
    plot(A, colorkey=FALSE, xlab='', ylab='', sub='') + 
    layer(panel.text(A@j+1L, A@i+1L, round(A@x, 2), 
                     col=gray(A@x>0.5))), ncol=3)

## ----aaa,include=FALSE---------------------------------------------------
proj <- inla.mesh.projector(mesh, dims=c(1.8,1)*200)
z <- inla.mesh.project(proj, field=c(0,0,0,0,0,0, 1, 0))
par(mar=c(0,0,0,0))
persp(proj$x, proj$y, z, xlab='x', ylab='y', theta=10, phi=80, col=gray(1), border=NA, shade=1/4, d=5)

