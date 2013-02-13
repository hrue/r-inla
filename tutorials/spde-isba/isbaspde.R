## This code accompanies a note in the ISBA Bulletin 2012, by
## Finn Lindgren, f.lindgren@bath.ac.uk
##
## Updated code and text available at
## http://people.bath.ac.uk/fl353/isba/
##
## 2013-01-14: Some minor code corrections. /FL
## 2013-01-16: Renamed covariate "offset" to "intercept" to avoid
##             confusion with the "offset()" formula feature. /FL
##
## Due to the complexity of builing the inla software, the R-INLA
## package is not available from CRAN, but can easily be installed
## using the following two lines of R code:
##
##   source("http://www.math.ntnu.no/inla/givemeINLA.R")
##   inla.upgrade(testing=TRUE)
##
## The second lines makes sure you get the latest development version,
## with the latest bug fixes and most up-to-date versions of the
## SPDE/GMRF model interface.

require(INLA)
require(lattice)

##############
m = 100
points = matrix(runif(m*2),m,2)
mesh = inla.mesh.create.helper(
  points=points,
  cutoff=0.05,
  offset=c(0.1,0.4),
  max.edge=c(0.05,0.5) )

## Plot the mesh and the data locations
plot(mesh)
points(points[,1],points[,2],pch=20,col="red")
## dev.copy2pdf(file="fig1.pdf")

##############
spde=inla.spde2.matern(mesh,alpha=2)

##############
sigma0 = 1   ## Standard deviation
range0 = 0.2 ## Spatial range
## Convert into tau and kappa:
kappa0 = sqrt(8)/range0
tau0 = 1/(sqrt(4*pi)*kappa0*sigma0)
## Construct spde object and set priors
## that will be used by inla().
spde=inla.spde2.matern(mesh,
  B.tau=cbind(log(tau0),1,0),
  B.kappa=cbind(log(kappa0),0,1),
  theta.prior.mean=c(0,0),
  theta.prior.prec=c(0.1,1))

##############
Q=inla.spde2.precision(
    spde,theta=c(0,0))
x=as.vector(inla.qsample(n=1,Q))
proj =
    inla.mesh.projector(mesh,
                        xlim=c(0,1),
                        ylim=c(0,1),
                        dims=c(200,200))
image(inla.mesh.project(
        proj,field=x))

##############
A=inla.spde.make.A(
  mesh,
  loc=points,
  index=rep(1:m,times=2),
  repl=rep(1:2,each=m) )

##############
Q=inla.spde.precision(
    spde,theta=c(0,0))
x=as.vector(inla.qsample(n=2,Q))
covariate = rnorm(m*2)
y = 5 + covariate*2 +
    as.vector(A %*% x) +
    rnorm(m*2)*0.01

##############
## Plot the first field replicate
cp = colorRampPalette(c("blue","cyan","yellow","red"))
at = pretty(c(-1,1)*max(abs(x[1:mesh$n])),101)
print(
levelplot(row.values=proj$x,
          column.values=proj$y,
          x=inla.mesh.project(
              proj,field=x[1:mesh$n]),
          xlim=c(0,1),
          ylim=c(0,1),
          col.regions=cp, at=at,
          aspect="iso",
          contour=FALSE, labels=FALSE, pretty=TRUE,
          xlab=NULL,ylab=NULL,scales=list(draw=FALSE))
    )
## dev.copy2pdf(file="fig2.pdf")

##############
mesh.index=inla.spde.make.index(
  name="field",
  n.mesh=mesh$n, n.repl=2)

##############
st.est=inla.stack(
  data=list(y=y),
  A=list(A,1),
  effects=list(
    c(mesh.index,list(intercept=1)),
    list(covar=covariate)),
  tag="est")

##############
st.pred=inla.stack(
  data=list(y=NA),
  A=list(1),
  effects=list(
    c(mesh.index,list(intercept=1))),
  tag="pred")

##############
stack = inla.stack(st.est,st.pred)

##############
formula =
  y ~ -1 + intercept + covar +
      f(field, model=spde,
        replicate=field.repl)
inla.result =
  inla(formula,
    data=inla.stack.data(stack,spde=spde),
    family="normal",
    control.predictor=
      list(A=inla.stack.A(stack),
           compute=TRUE))

##############
result = inla.spde2.result(
  inla.result, "field", spde)
plot(result$marginals.range.nominal[[1]])

##############
index=inla.stack.index(
  stack,"pred")$data
image(inla.mesh.project(proj,
  inla.result$summary.linear.predictor$mean[
    index[mesh.index$field.repl==1]]))
image(inla.mesh.project(proj,
  inla.result$summary.linear.predictor$sd[
    index[mesh.index$field.repl==1]]))

##############
## Posterior mean
print(
levelplot(row.values=proj$x,
          column.values=proj$y,
          x=inla.mesh.project(proj,
          inla.result$summary.linear.predictor$mean[index[mesh.index$field.repl==1]]),
          xlim=c(0,1),
          ylim=c(0,1),
          col.regions=cp, at=at+5,
          aspect="iso",
          contour=FALSE, labels=FALSE, pretty=TRUE,
          xlab=NULL,ylab=NULL,scales=list(draw=FALSE))
    )
## dev.copy2pdf(file="fig3.pdf")

## Posterior std.dev.
print(
levelplot(row.values=proj$x,
          column.values=proj$y,
          x=inla.mesh.project(proj,
          inla.result$summary.linear.predictor$sd[index[mesh.index$field.repl==1]]),
          xlim=c(0,1),
          ylim=c(0,1),
          col.regions=cp,
          aspect="iso",
          contour=FALSE, labels=FALSE, pretty=TRUE,
          xlab=NULL,ylab=NULL,scales=list(draw=FALSE))
    )
## dev.copy2pdf(file="fig4.pdf")
