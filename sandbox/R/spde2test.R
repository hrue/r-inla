remove(list=ls())
graphics.off()
require(rgl)
rgl.quit()
require(INLA)
inla.my.update()

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

spde1 = inla.spde.create(mesh)
spde1$f$hyper.default =
    list(theta1 = list(initial=0, param=c(0,prior.prec)),
         theta2 = list(initial=0, param=c(0,prior.prec)),
         theta3 = list(initial=0, param=c(0,prior.prec))
         )

##require(debug)
##mtrace(inla.spde.create.generic)
spde2 =
    inla.spde.generic2(M0=spde1$internal$c0,
                       M1=spde1$internal$g1,
                       M2=spde1$internal$g2,
                       B0=matrix(c(0,1,0),1,3),
                       B1=matrix(c(0,0,1),1,3),
                       B2=1,
                       theta.mu = c(0,0),
                       theta.Q = diag(nrow=2)*prior.prec,
                       transform="identity")

spde3 = inla.spde.matern(mesh)

field = inla.spde.query(spde1, sample=list(tau=tau, kappa2=kappa2))$sample
Y = 1+rep(field, n.obs.repl)+rnorm(mesh$n*n.obs.repl)/sqrt(prec)

data1 = list(Y=Y, field=rep(1:mesh$n, n.obs.repl), model=spde1)
formula1 = (as.formula("Y ~ 1 + f(field, model=spde1)"))
data2 = list(Y=Y, field=rep(1:mesh$n, n.obs.repl), model=spde2)
formula2 = (as.formula("Y ~ 1 + f(field, model=spde2)"))
data3 = list(Y=Y, field=rep(1:mesh$n, n.obs.repl), model=spde3)
formula3 = (as.formula("Y ~ 1 + f(field, model=spde3)"))

verbose = TRUE
result1 =
    inla(formula=formula1,
         family="gaussian",
         data=data1,
         control.predictor=list(compute=TRUE),
         control.compute=list(cpo=TRUE),
         verbose=verbose)
result2 =
    inla(formula=formula2,
         family="gaussian",
         data=data2,
         control.predictor=list(compute=TRUE),
         control.compute=list(cpo=TRUE),
         verbose=verbose)
result3 =
    inla(formula=formula3,
         family="gaussian",
         data=data3,
         control.predictor=list(compute=TRUE),
         control.compute=list(cpo=TRUE),
         verbose=verbose)

result =
    cbind(theta,
          result1$summary.hyperpar[,c("0.025quant","0.975quant")],
          result2$summary.hyperpar[,c("0.025quant","0.975quant")],
          rbind(NA, result2$summary.spde2.blc$field[,c("0.025quant","0.975quant")]))
colnames(result) <- c("Truth", rep("SPDE1", 2), rep("SPDE2", 2), rep("SPDE2(BLC)", 2))
rownames(result) <- c("Prec.", "theta1", "theta2")
print(result)

require(fields)
plot(old.mesh.class(mesh), field,
     color.palette=tim.colors,
     draw.edges=FALSE, draw.vertices=FALSE)
plot(old.mesh.class(mesh), result2$summary.random$field[,"mean"],
     color.palette=tim.colors,
     draw.edges=FALSE, draw.vertices=FALSE)
plot(old.mesh.class(mesh), result2$summary.random$field[,"sd"],
     color.palette=tim.colors,
     draw.edges=FALSE, draw.vertices=FALSE)


ks.plot = function (x, y, ...) {
    test = ks.test(x, y, ...)
    n = length(x)
    Fn = ((1:n)-0.5)/n
    F = y(sort(x))
    empirical.diff = (Fn-F)*sqrt(n)
    T = max(abs(empirical.diff))
    ylim = c(-1,1)*max(1,T)
    plot(Fn, empirical.diff, type='l',
         ylim=ylim,
         main=paste("K-S-test, p-value = ",test$p.value),
         ylab="(Fn-F) sqrt(n)",
         xlab="Quantile"
         )
    lines(Fn, 2*sqrt(Fn*(1-Fn)), type='l')
    lines(Fn, -2*sqrt(Fn*(1-Fn)), type='l')
    lines(c(0,1), c(1,1)*T, type='l')
    lines(c(0,1), -c(1,1)*T, type='l')
    invisible(test)
}

dev.new()
ks.plot(result1$pit, punif)
dev.new()
ks.plot(result2$pit, punif)
dev.new()
ks.plot(result3$pit, punif)



if (FALSE) {
## Option 1:
spde3 = inla.spde.create(mesh, model="matern", param=list(alpha=2, basis.T=...))
Q.epsilon = inla.spde.query(spde3, precision=list(tau=1, kappa2=0.5))$precision
spde4 = inla.spde.create(mesh, model="heatequation",
   param=list(Q.epsilon=Q.epsilon))
## Option 1.5:
spde3 = inla.spde.create(mesh, model="matern", alpha=2, basis.T=...)
Q.epsilon = inla.spde.query(spde3, precision=TRUE, tau=1, kappa2=0.5)$precision
spde4 = inla.spde.create(mesh, model="heatequation", Q.epsilon=Q.epsilon)
## Option 2:
spde3 = inla.spde.create.matern(mesh, alpha=2)
Q.epsilon = inla.spde.precision(spde3, tau=1, kappa2=0.5)
spde4 = inla.spde.create.heatequation(mesh, Q.epsilon=Q.epsilon)
}
