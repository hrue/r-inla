n = 100
idx = 1:n
y = sin((idx/n)^3 * 4*pi) + rnorm(n, sd=0.2)

inla.setOption(scale.model.default = FALSE)
log.prec = 2
dev.new()
plot(idx, y, main="default, scale.model=F")
r1 = inla(y ~ 1 + f(idx, model="rw1", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
r2 = inla(y ~ 1 + f(idx, model="rw2", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
lines(idx,  r1$summary.linear.predictor$mean, lwd=1)
lines(idx,  r2$summary.linear.predictor$mean, lwd=2)

inla.setOption(scale.model.default = FALSE)
log.prec = 2
dev.new()
idx = (1:n)/100
plot(idx, y, main="default, scale.model=F, idx /= 100")
r1 = inla(y ~ 1 + f(idx, model="rw1", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
r2 = inla(y ~ 1 + f(idx, model="rw2", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
lines(idx,  r1$summary.linear.predictor$mean, lwd=1)
lines(idx,  r2$summary.linear.predictor$mean, lwd=2)

inla.setOption(scale.model.default = TRUE)
dev.new()
idx = 1:n
plot(idx, y, main="default, scale.model=T")
r1 = inla(y ~ 1 + f(idx, model="rw1", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
r2 = inla(y ~ 1 + f(idx, model="rw2", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
lines(idx,  r1$summary.linear.predictor$mean, lwd=1)
lines(idx,  r2$summary.linear.predictor$mean, lwd=2)

inla.setOption(scale.model.default = TRUE)
dev.new()
idx = (1:n)/100
plot(idx, y, main="default, scale.model=T,idx /= 100")
r1 = inla(y ~ 1 + f(idx, model="rw1", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
r2 = inla(y ~ 1 + f(idx, model="rw2", initial=log.prec, fixed=TRUE),
        data = data.frame(y, idx),
        control.predictor = list(compute=TRUE))
lines(idx,  r1$summary.linear.predictor$mean, lwd=1)
lines(idx,  r2$summary.linear.predictor$mean, lwd=2)
