## CGD: Chronic granulomatous disease, see Manda and Meyer (2005)

use.hyperpar = TRUE
write.new.figures = TRUE
plot.figures = TRUE
new.windows = TRUE
bty = "l"
lwd = 2L

library(survival)
data(cgd)

## Hospital region: US versus Europe
cgd$region[cgd$hos.cat == "US:NIH"] <- 0
cgd$region[cgd$hos.cat == "US:other"] <- 0
cgd$region[cgd$hos.cat == "Europe:Amsterdam"] <- 1
cgd$region[cgd$hos.cat == "Europe:other"] <- 1

## Pattern of inheritance
cgd$inherit2[cgd$inherit == "X-linked"] <- 1
cgd$inherit2[cgd$inherit == "autosomal"] <- 0

cgd$time = cgd$tstop-cgd$tstart

## we find the time for the first event for each patient. set it to -1
## if no event. we also need to define the variable 'deterministic' to
## be completed later
cgd$deterministic = 0
cgd$first.time = -1L
for (id in unique(cgd$id)) {
    idx = which(cgd$id == id)
    if (any(cgd[idx, "status"] == 1)) {
        ## there is an infection for person id
        idxx = which(cgd[idx, "status"] == 1)
        first.time = min(cgd[idxx, "tstop"])
        cgd[idx, "first.time"] = first.time
    }
}

## Number of intervals for piecewise hazard function
nc = 25

## For the rho1 prior
upper.rho = 1/2
alpha.rho = 0.75 ## lower limit = 1/2

## For the spline priors (gives stdev = 0.31*u.)
u.frailty = 0.30/0.31
a.frailty = 0.01
u.bh = 0.15/0.31
a.bh = 0.01

## scale numerical covariates for the fixed effects
cgd$age = scale(cgd$age)
cgd$height = scale(cgd$height)
cgd$weight = scale(cgd$weight)

p = inla.coxph(inla.surv(time, status) ~ 1 +
        treat + inherit2 + age + height + weight + propylac + sex + region +
        steroids + deterministic + 
        f(baseline.hazard.idx, model = "ar1", replicate = id, values = 1:nc, 
          hyper = list(
                  prec = list(
                          prior = "pc.prec",
                          param = c(u.frailty, a.frailty)),
                  rho = list(
                          prior = "pc.rho1",
                          param = c(upper.rho, alpha.rho)))),
        data = cgd,
        control.hazard = list(
                model = "rw1",
                n.intervals = nc,
                scale.model = T,
                constr = TRUE, 
                hyper = list(
                        prec = list(
                                prior = "pc.prec",
                                param = c(u.bh, a.bh)))))

## now that the data.frame is expanded, we need to compute the
## 'deterministic' variable, which is the positive time from the first
## event. the time for the first event is 'first.time'
for (id in unique(p$data$id)) {
    idx = which(p$data$id == id)
    if (any(p$data[idx, "first.time"] >= 0)) {
        ## there is an event, and all $first.time should be
        ## equal. internal check:
        stopifnot(all(p$data[idx, "first.time"] == p$data[idx[1], "first.time"]))

        ## we only take the positive time, meaning that the time
        ## before the first.time is taken as 0
        p$data[idx, "deterministic"] = pmax(0,
                      p$data[idx, "tstart"] +
                      p$data[idx, "baseline.hazard.time"] -
                      p$data[idx, "first.time"])
    }
}
p$data$deterministic = p$data$deterministic / 365  ## this is just a
                                                   ## rescaling...
result = inla(p$formula, 
        family = p$family, 
        data=c(as.list(p$data), p$data.list), 
        control.fixed = list(
                prec.intercept = 0,
                prec = 0.001),
        E = p$E,
        verbose = TRUE)

if (use.hyperpar) 
    result = inla.hyperpar(result, diff.logdens = 10, verbose=TRUE)

## rho
if (plot.figures) {
    if (new.windows)
        inla.dev.new()
    m = inla.smarginal(
        inla.tmarginal(function(x) 2 * 1/(1+exp(-x)) -1,
                       inla.smarginal(
                           result$internal.marginals.hyperpar[["Rho_intern for baseline.hazard.idx"]],
                           extrapolate = 1,  factor = 100)),
        extrapolate = 1,  factor = 100)
    idx = which(m$x < 1.0 & m$x > 0)
    m = list(x = m$x[idx], y = m$y[idx])
    plot(m, 
         type="l",  lwd=lwd, bty = bty,
         log = "y", 
         xlab = "Lag-one correlation", ylab = "Density")

    lines(m$x, inla.pc.dcor1(m$x, upper.rho, alpha.rho), lwd=lwd, lty=2)
    if (write.new.figures) {
        dev.print(postscript,  file = "coxph-rho.ps")
    }

    ## same plot not it log-scale
    if (new.windows)
        inla.dev.new()
    plot(m, 
         type="l", lwd=lwd, bty = bty, ylim = c(0,10),
         xlab = "Lag-one correlation", ylab = "Density")
    lines(m$x, inla.pc.dcor1(m$x, upper.rho, alpha.rho), lwd=lwd, lty=2)
    if (write.new.figures) {
        dev.print(postscript,  file = "coxph-rho-2.ps")
    }
}

## prec for baseline hazard
m = inla.tmarginal(exp,
        inla.smarginal(result$internal.marginals.hyperpar[["Log precision for baseline.hazard"]],
                       extrapolate = 3, factor = 1000L))
idx = which(m$x < 2000)
m = list(x = c(0, m$x),  y = c(.Machine$double.eps, m$y))
m = inla.smarginal(list(x = m$x[idx],  y = m$y[idx]),  factor = 1000L, extrapolate = 1)
idx = which(m$x <= 0)
m = list(x = m$x[-idx],  y = m$y[-idx])

if (plot.figures) {
    if (new.windows)
        inla.dev.new()
    plot(m,
         type = "l",
         xlim = c(0, 300), 
         ylim = c(0, 0.008), 
         lwd = lwd,  bty = bty,
         xlab = "Precision for the log.baseline.hazard",
         ylab = "Density")
    prec = seq(0.001, 400, len=10000)
    lines(prec, inla.pc.dprec(prec, u.bh, a.bh, log=FALSE), lty=2, lwd=2)
    if (write.new.figures) 
        dev.print(postscript,  file = "coxph-bh-prec.ps")

    
    if (new.windows)
        inla.dev.new()
    xx = p$data.list$baseline.hazard.values
    plot(NA, NA, xlim = range(xx), ylim = c(-0.6, 0.8),
         xlab = "Relative time",
         ylab = "Log baseline.hazard",
         bty = bty)
    lines(xx, result$summary.random$baseline.hazard$mean,
          type="s", lwd = lwd)
    lines(xx, result$summary.random$baseline.hazard$"0.025quant", 
          type="s", lwd = lwd,  lty = 2)
    lines(xx, result$summary.random$baseline.hazard$"0.975quant", 
          type="s", lwd = lwd,  lty = 3)
    lines(xx, result$summary.random$baseline.hazard$"0.5quant", 
          type="s", lwd = lwd,  lty = 4)
    abline(a=0, b=0)
    if (write.new.figures) 
        dev.print(postscript,  file = "coxph-bh.ps")

    print(paste("posterior sigma*",  exp(mean(log(result$summary.random$baseline.hazard$sd)))))
    print(paste("prior sigma*", 0.3 * u.bh))
}

if (write.new.figures) {
    system("which psfix && psfix coxph*.ps")
}
