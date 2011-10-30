remove(list=ls())
graphics.off()
require(rgl)
rgl.quit()
require(INLA)
##inla.my.update()

n.obs.repl = 2
prec = 10^2
sigma2 = 1
kappa2 = 8/0.2^2
tau = 1/sqrt(4*pi*kappa2*sigma2)
theta = c(prec, log(tau), log(kappa2))
prior.prec = 0.5

##mesh =
##    inla.mesh.create(boundary=matrix(c(0,0,1,0,1,1,0,1),4,2,byrow=TRUE),
##                     refine=list(max.edge=0.05),
##                     plot.delay=0)
mesh =
    inla.mesh.create(lattice=inla.mesh.lattice(dims=c(40,40)),
                     extend=FALSE,
                     plot.delay=NULL)

spde1 = inla.spde1.create(mesh)
spde1$f$hyper.default =
    list(theta1 = list(initial=0, param=c(0,prior.prec)),
         theta2 = list(initial=0, param=c(0,prior.prec)),
         theta3 = list(initial=0, param=c(0,prior.prec))
         )

##require(debug)
##mtrace(inla.spde.create.generic)
spde2a =
    inla.spde2.generic(M0=spde1$internal$c0,
                       M1=spde1$internal$g1,
                       M2=spde1$internal$g2,
                       B0=matrix(c(0,1,0),1,3),
                       B1=matrix(c(0,0,1),1,3),
                       B2=1,
                       theta.mu = c(0,0),
                       theta.Q = diag(nrow=2)*prior.prec,
                       transform="identity")

spde2b = inla.spde2.matern(mesh)

field = inla.spde1.query(spde1, sample=list(tau=tau, kappa2=kappa2))$sample
Y = 1+rep(field, n.obs.repl)+rnorm(mesh$n*n.obs.repl)/sqrt(prec)

data1 = list(Y=Y, field=rep(1:mesh$n, n.obs.repl), spde=spde1)
formula1 = (as.formula("Y ~ 1 + f(field, model=spde)"))
data2 = list(Y=Y, field=rep(1:mesh$n, n.obs.repl), spde=spde2a)
formula2 = (as.formula("Y ~ 1 + f(field, model=spde)"))
data3 = list(Y=Y, field=rep(1:mesh$n, n.obs.repl), spde=spde2b)
formula3 = (as.formula("Y ~ 1 + f(field, model=spde)"))

verbose = FALSE
result1 =
    inla(formula=formula1,
         family="gaussian",
         data=data1,
         control.predictor=list(compute=TRUE),
         control.compute=list(cpo=TRUE),
         verbose=verbose)
result2a =
    inla(formula=formula2,
         family="gaussian",
         data=data2,
         control.predictor=list(compute=TRUE),
         control.compute=list(cpo=TRUE),
         verbose=verbose)
result2b =
    inla(formula=formula3,
         family="gaussian",
         data=data3,
         control.predictor=list(compute=TRUE),
         control.compute=list(cpo=TRUE),
         verbose=verbose)

BLC.names = rownames(spde2b$internal$param.generic$BLC)
rownames(result2b$summary.spde2.blc$field) <- BLC.names
names(result2b$marginals.spde2.blc$field) <-  BLC.names

result =
    cbind(theta,
          result1$summary.hyperpar[,c("0.025quant","0.975quant")],
          result2a$summary.hyperpar[,c("0.025quant","0.975quant")],
          rbind(NA, result2a$summary.spde2.blc$field[1:2,c("0.025quant","0.975quant")]),
          result2b$summary.hyperpar[,c("0.025quant","0.975quant")],
          rbind(NA, result2b$summary.spde2.blc$field[1:2,c("0.025quant","0.975quant")]))
colnames(result) <- c("Truth", "SPDE", rep("SPDE2A", 2), rep("SPDE2A(BLC)", 2), rep("SPDE2B", 2), rep("SPDE2B(BLC)", 2))
rownames(result) <- c("Prec.", "theta1", "theta2")
print(result[2:3,])

if (FALSE) {
require(fields)
plot(old.mesh.class(mesh), field,
     color.palette=tim.colors,
     draw.edges=FALSE, draw.vertices=FALSE)
plot(old.mesh.class(mesh), result2a$summary.random$field[,"mean"],
     color.palette=tim.colors,
     draw.edges=FALSE, draw.vertices=FALSE)
plot(old.mesh.class(mesh), result2a$summary.random$field[,"sd"],
     color.palette=tim.colors,
     draw.edges=FALSE, draw.vertices=FALSE)
}

##dev.new()
##inla.ks.plot(result1$pit, punif)
dev.new()
inla.ks.plot(result2a$pit, punif)
dev.new()
inla.ks.plot(result2b$pit, punif)

dev.new()
plot(result2a$marginals.hyperpar[[2]],
     main="SPDE2A, theta.1: hyperpar (circles), spde2.blc (solid)")
lines(result2a$marginals.spde2.blc$field[[1]])
##dev.new()
##plot(result2a$marginals.hyperpar[[3]],
##     main="SPDE2A, theta.2: hyperpar (circles), spde2.blc (solid)")
##lines(result2a$marginals.spde2.blc$field[[2]])
dev.new()
plot(result2b$marginals.hyperpar[[2]],
     main="SPDE2B, theta.1: hyperpar (circles), spde2.blc (solid)")
lines(result2b$marginals.spde2.blc$field[[1]])
dev.new()
plot(result2b$marginals.hyperpar[[3]],
     main="SPDE2B, theta.2: hyperpar (circles), spde2.blc (solid)")
lines(result2b$marginals.spde2.blc$field[[2]])


if (FALSE) {
## Option 1:
spde2b = inla.spde.create(mesh, model="matern", param=list(alpha=2, basis.T=...))
Q.epsilon = inla.spde.query(spde2b, precision=list(tau=1, kappa2=0.5))$precision
spde4 = inla.spde.create(mesh, model="heatequation",
   param=list(Q.epsilon=Q.epsilon))
## Option 1.5:
spde2b = inla.spde.create(mesh, model="matern", alpha=2, basis.T=...)
Q.epsilon = inla.spde.query(spde2b, precision=TRUE, tau=1, kappa2=0.5)$precision
spde4 = inla.spde.create(mesh, model="heatequation", Q.epsilon=Q.epsilon)
## Option 2:
spde2b = inla.spde.matern(mesh, alpha=2)
Q.epsilon = inla.spde.precision(spde2b, tau=1, kappa2=0.5)
spde4 = inla.spde.heatequation(mesh, Q.epsilon=Q.epsilon)
}
