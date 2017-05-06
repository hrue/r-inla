## ----settings,include=FALSE,message=FALSE, warning=FALSE-----------------
library(knitr)
opts_chunk$set(
fig.path='figs/lowst',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ") 
library(splancs) 
library(INLA) 
library(fields) 
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)

## ----smesh---------------------------------------------------------------
data(PRprec)
bound <- inla.nonconvex.hull(as.matrix(PRprec[,1:2]), .2, .2, resol=50)
mesh.s <- inla.mesh.2d(bound=bound, max.edge=c(1,2),
                       offset=c(1e-5,0.7), cutoff=0.5)
spde.s <- inla.spde2.matern(mesh.s)

## ----dat-----------------------------------------------------------------
dim(PRprec)
PRprec[1:2, 1:7]

## ----datagg--------------------------------------------------------------
table(table(id5 <- 0:364%/%5 + 1))
n5 <- t(apply(!is.na(PRprec[,3+1:365]), 1, tapply, id5, sum))
y5 <- t(apply(PRprec[,3+1:365]>0.1, 1, tapply, id5, sum, na.rm=TRUE))
k <- ncol(n5);       table(as.vector(n5))

## ----na5-----------------------------------------------------------------
y5[n5==0] <- NA;     n5[n5==0] <- 5

## ----tmesh---------------------------------------------------------------
bt <- 6;   gtime <- seq(1+bt, k, length=round(k/bt))-bt/2
mesh.t <- inla.mesh.1d(gtime, degree=1)
table(igr <- apply(abs(outer(mesh.t$loc, 1:k, '-')), 2, which.min))

## ----nk------------------------------------------------------------------
spde.s$n.spde*mesh.t$n

## ----repcoo--------------------------------------------------------------
n <- nrow(PRprec)
st.sloc <- cbind(rep(PRprec[,1], k), rep(PRprec[,2], k))

## ----Ast-----------------------------------------------------------------
Ast <- inla.spde.make.A(mesh=mesh.s, loc=st.sloc, 
                        group.mesh=mesh.t, group=rep(1:k, each=n))

## ----stk-----------------------------------------------------------------
idx.st <- inla.spde.make.index('i', n.spde=spde.s$n.spde,
                               n.group=mesh.t$n)
dat <- inla.stack(data=list(yy=as.vector(y5), nn=as.vector(n5)), 
                  A=list(Ast, 1), 
                  effects=list(idx.st, 
                      data.frame(mu0=1, 
                                 altitude=rep(PRprec$Alt/1e3, k))))

## ----form----------------------------------------------------------------
form <- yy ~ 0 + mu0 + altitude + 
    f(i, model=spde.s, group=i.group,
      control.group=list(model='ar1'))

## ----res,results='hide'--------------------------------------------------
result <- inla(form, 'binomial', data=inla.stack.data(dat),
               Ntrials=inla.stack.data(dat)$nn, 
               control.predictor=list(A=inla.stack.A(dat)),
               control.mode=list(theta=c(-0.48, -0.9, 2.52)), ###restart=TRUE), 
               control.inla=list(strategy='gaussian', int.strategy='eb'))

## ----grid----------------------------------------------------------------
data(PRborder)
r0 <- diff(range(PRborder[,1]))/diff(range(PRborder[,2]))
prj <- inla.mesh.projector(mesh.s, xlim=range(PRborder[,1]),
                           ylim=range(PRborder[,2]), dims=c(100*r0, 100))
in.pr <- inout(prj$lattice$loc, PRborder)

## ----proj----------------------------------------------------------------
mu.spat <- lapply(1:mesh.t$n, function(j) {
  r <- inla.mesh.project(prj, field=result$summary.ran$i$mean[
                         1:spde.s$n.spde + (j-1)*spde.s$n.spde])
  r[!in.pr] <- NA;   return(r)})

## ----plt,eval=F----------------------------------------------------------
## par(mfrow=c(4,3), mar=c(0,0,0,0))
## zlm <- range(unlist(mu.spat), na.rm=TRUE)
## for (j in 1:mesh.t$n) {
##     image.plot(x=prj$x, y=prj$y, z=mu.spat[[j]], asp=1, axes=FALSE, zlim=zlm)
##     lines(PRborder)
##     points(PRprec[, 1:2],
##            cex=rowSums(y5[, j==igr], na.rm=TRUE)/rowSums(n5[,j==igr]))
## }

## ----lowstres,echo=FALSE,fig.width=10,fig.height=8.7---------------------
par(mfrow=c(4,3), mar=c(0,0,0,0))
zlm <- range(unlist(mu.spat), na.rm=TRUE)
for (j in 1:mesh.t$n) {
    image.plot(x=prj$x, y=prj$y, z=mu.spat[[j]], asp=1, axes=FALSE, zlim=zlm)
    lines(PRborder)
    points(PRprec[, 1:2], 
           cex=rowSums(y5[, j==igr], na.rm=TRUE)/rowSums(n5[,j==igr]))
}

