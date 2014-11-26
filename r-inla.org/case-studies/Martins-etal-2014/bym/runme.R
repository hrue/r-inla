library(INLA)
write.new.figures = TRUE

inla.dev.new.hook = function()
{
    cex.lab = 1.4
    cex.axis = 1.4
    par(cex.lab=cex.lab, cex.axis=cex.axis)
}
inla.dev.new()

g = "sardinia.graph"
sardinia = read.table("sardinia.dat", col.names=c("y", "E", "SMR"))
n = dim(sardinia)[1]
sardinia$idx = 1:n
Q = INLA:::inla.pc.bym.Q(g)
Q = INLA:::inla.scale.model(Q,  constr=list(A=matrix(1, 1, n), e=0))
u = 0.2/0.31
alpha = 0.01
phi.u = 0.5
phi.alpha = 2/3 ## prob(phi < phi.u) = phi.alpha

if (FALSE) {
    ## make figures of the priors
    for (g in c("germany.graph", "sardinia.graph")) {
        log.prior = INLA:::inla.pc.bym.phi(Q=Q, rankdef=1, u=phi.u, alpha = phi.alpha)
        phis = 1/(1+exp(-seq(-10, 10, len=10000)))

        inla.dev.new()
        plot(phis, exp(log.prior(phis)), type="l",
             lwd = 2,  bty = "l",
             xlab = expression(phi), 
             ylab = "Density")
        if (length(grep("^sar", g)) > 0) {
            fnm = "sardinia-prior.ps"
        } else {
            fnm = "germany-prior.ps"
        }
        if (write.new.figures)
            dev.print(postscript, file=fnm)
    }
    if (write.new.figures)
        system("which psfix && psfix *-prior.ps")
}

phis = 1/(1+exp(-seq(-15, 15,  len = 10000)))
phi.prior = INLA:::inla.pc.bym.phi(Q=Q, u= phi.u, alpha = phi.alpha)

if (FALSE) {
    cat("\n\nSIMULATE NEW DATA\n\n\n")
    s = 0
    if (TRUE) {
        ## simulate perfect iid-data
        sardinia$y = rpois(n, lambda = sardinia$E * exp(rnorm(n, sd = s)))
    } else {
        ## simulate perfect spatial data with no unstructured term
        QQ = Q
        diag(QQ) = diag(QQ) + 1e-6
        x = inla.qsample(1, Q=QQ, constr = list(A = matrix(1, 1, n), e=0))
        sardinia$y = rpois(n, lambda = sardinia$E * exp(s*x))
    }
}

formula = y ~ 1 + f(idx,
        model = "bym2",
        graph=g,
        scale.model = TRUE, 
        constr = TRUE, 
        hyper=list(
                phi = list(
                        prior = "pc",
                        param = c(phi.u, phi.alpha),
                        initial = -3), 
                prec = list(
                        prior = "pc.prec",
                        param = c(u, alpha),
                        initial = 5)))

r = inla(formula, data = sardinia, family = "poisson", E=E,
        control.predictor = list(compute=TRUE))
r = r.sardinia = inla.hyperpar(r, dz = 0.2, diff.logdens=20)

summary(r)

## the maps
source("sardinia-draw.R")
inla.dev.new()
sardinia.map(r$summary.fitted.values$mean)
title("Relative risk")
if (write.new.figures)
    dev.print(postscript, file = "sardinia-relative-risk.ps")

## prec
m.p = inla.smarginal(r$internal.marginals.hyperpar[[
        "Log precision for idx"]], factor = 100, extrapolate = 0.5)
mm.p = inla.tmarginal(function(x) exp(x), m.p)
inla.dev.new()
plot(mm.p, type="l", lwd = 2, bty="l", 
     xlim = c(0, 300),
     xlab = "Precision",
     ylab = "Density")
prec = seq(0.001, 3000, len=10000)
lines(prec, inla.pc.dprec(prec, u, alpha, log=FALSE), lty=2, lwd=2)
if (write.new.figures)
    dev.print(postscript, file = "sardinia-precision.ps")

nn = 100000
x = rnorm(nn, s = sqrt(1/exp(inla.rmarginal(nn, m.p))))
s = sd(x)
print(paste("marginal stdev", s))

## phi
m.r = inla.smarginal(r$internal.marginals.hyperpar[[
        "Logit phi for idx"]], factor = 100, extrapolate = 0.5)
mm.r = inla.tmarginal(function(x)1/(1+exp(-x)), m.r)
inla.dev.new()
plot(mm.r, type="l", lwd = 2, bty="l",
     xlim = c(0, 1),
     ylim = c(0, 10),
     xlab = expression(phi), 
     ylab = "Density")
lines(phis, exp(phi.prior(phis)), lwd = 2,  lty=2)
if (write.new.figures)
    dev.print(postscript, file = "sardinia-phi.ps")

if (write.new.figures)
    system("which psfix && psfix sardinia-*.ps")


#####################################################################
## the Germany example

data(Germany)
g = "germany.graph"
n = dim(Germany)[1]
Q = INLA:::inla.pc.bym.Q(g)
Q = INLA:::inla.scale.model(Q, constr=list(A=matrix(1, 1, n), e=0))

u = 0.2/0.31
alpha = 0.01
phi.u = 0.5
phi.alpha = 2/3 ## prob(phi < phi.u) = phi.alpha

phi.prior = INLA:::inla.pc.bym.phi(Q=Q, u= phi.u, alpha = phi.alpha)
values = (-9):110
values = 1:100
formula = Y ~ 1 + f(region,
        model = "bym2",
        graph=g,
        scale.model = TRUE, 
        constr = TRUE, 
        hyper=list(
                phi = list(
                        prior = "pc",
                        param = c(phi.u, phi.alpha), 
                        initial = -3), 
                prec = list(
                        prior = "pc.prec",
                        param = c(u, alpha),
                        initial = 5))) +
    f(x, model="rw2", values = values, 
      scale.model=TRUE, 
      hyper = list(
              prec = list(
                      prior = "pc.prec",
                      param = c(u, alpha))))
                      

r = inla(formula, data = Germany, family = "poisson", E=E,
        control.predictor = list(compute = TRUE))
r = r.germany = inla.hyperpar(r, dz = .5, diff.logdens=15)
summary(r)

if (exists("germany.map")) rm("germany.map")
source("germany-draw.R")
germany.map(r$summary.fitted.values$mean)
title("Relative risk")
if (write.new.figures)
    dev.print(postscript, file = "germany-relative-risk.ps")

## prec
m.p = inla.smarginal(r$internal.marginals.hyperpar[[
        "Log precision for region"]],
        factor = 100, extrapolate = .5)
mm.p = inla.tmarginal(function(x) exp(x), m.p)
mm.p = inla.smarginal(mm.p, extrapolate = .25)
inla.dev.new()
plot(mm.p, type="l", lwd = 2, bty="l",
     xlim = c(0, 80),
     xlab = "Precision",
     ylab = "Density")
prec = seq(0.001, 200, len=10000)
lines(prec, inla.pc.dprec(prec, u, alpha, log=FALSE), lty=2, lwd=2)
if (write.new.figures)
    dev.print(postscript, file = "germany-precision.ps")

nn = 100000
x = rnorm(nn, s = sqrt(1/exp(inla.rmarginal(nn, m.p))))
s = sd(x)
print(paste("marginal stdev", s))

## phi
m.r = inla.smarginal(r$internal.marginals.hyperpar[["Logit phi for region"]],
        factor = 100, extrapolate = 0.5)
mm.r = inla.tmarginal(function(x)1/(1+exp(-x)), m.r)
mm.r = inla.smarginal(mm.r, extrapolate = .6)
inla.dev.new()
plot(mm.r, type="l", lwd = 2, bty="l",
     xlim = c(0, 1),
     ylim = c(0, 10),
     xlab = expression(phi), 
     ylab = "Density")
lines(phis, exp(phi.prior(phis)), lwd = 2,  lty=2)
if (write.new.figures)
    dev.print(postscript, file = "germany-phi.ps")

## covariate
m.x = r$summary.random$x
inla.dev.new()
plot(m.x[, "ID"], m.x[, "mean"], type="l",  lwd = 3, bty="l",
     xlim = c(min(values)-1, max(values)+1),
     ylim = c(-0.6, 0.6),
     xlab = "Covariate",
     ylab = "f(covariate)")
lines(m.x[, "ID"], m.x[, "0.025quant"], lwd = 2, lty = 3)
lines(m.x[, "ID"], m.x[, "0.975quant"], lwd = 2, lty = 3)

yy = m.x[, "mean"]
xx = m.x[, "ID"]
reg = lm(yy ~ 1 + xx)
lines(xx, coef(reg)[1] + coef(reg)[2]*xx, 
      lwd=1, lty=2, col = "blue")

if (write.new.figures) {
    dev.print(postscript, file = "germany-covariate.ps")
    system("which psfix && psfix germany-*.ps")
}
