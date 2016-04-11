## set variables run.inla=TRUE/FALSE, if TRUE start inla() if FALSE run the R with the prior

if (!exists("run.inla")) run.inla = FALSE

source("aux-core.R")
source("pcp-tools.R")
source("ranf.R")
library(numDeriv)
library(mvtnorm)
library(INLA)
inla.setOption(scale.model.default=TRUE)

make.constr = function(values)
{
    ## make the constraints for RW2 with given (integer) values
    n = length(values)
    A = matrix(NA, 2, n)
    return (list(A=rbind(1, 1:n - (n+1)/2), e = rep(0, 2)))
}

## the big logreg example
library("numDeriv")
data("heart", package = "catdata")
heart = as.data.frame(heart)
##heart = heart[1:200, ]
n = dim(heart)[1]

## Adopt the variable names as in the Wood-Kohn paper, and make values go from 1...n
heart$CR = heart$CR2 = round(heart$ldl*10) - min(round(heart$ldl*10)) + 1L
heart$BP = heart$BP2 = round(heart$sbp) - min(round(heart$sbp)) + 1L
heart$Age = heart$Age2 = round(heart$age) - min(round(heart$age)) + 1L
CR.values = seq(min(heart$CR), max(heart$CR), by = 1L)
BP.values = seq(min(heart$BP), max(heart$BP), by = 1L)
Age.values = seq(min(heart$Age), max(heart$Age), by = 1L)

heart$BP.idx = rep(1, n)
heart$CR.idx = rep(1, n)
heart$Age.idx = rep(1, n)

pc.prec.param = c(2, 0.01)
diag.value = 1e-6



model = as.environment(
    list(
        n.terms = 3,
        delta.add.diag = 1e-6, 
        n = dim(heart)[1], 
        X = list(BP=heart$BP, CR=heart$CR, Age=heart$Age), 
        S.fixed = list(), 
        S.random = list(), 
        param.case1 = list(U=0.5, alpha=0.75), 
        param.case2 = list(lambda = 0.3), 
        param.prec = list(U=4.84, alpha = 0.01), 
        param.intercept = list(variance = 100),

        ## hold case priors. these do not not need to be recomputed.
        priors.case1 = list(),

        ## hold case2 prior. the case2 prior depends on hyperparameters and must be recomputed
        ## all the time. Maybe its faily constant?
        priors.case2 = list(
            old = c(), 
            new = c()), 

        ## weights between the contributions from the different covariates. the mapping is given
        ## by the w2phi() and phi2w() functions.
        weights.intern = list(
            old = c(),
            new = c()), 
        weights = list(
            old = c(),
            new = c()), 

        ## variance contribution from the spline. phi=0 is the linear effect only, and phi=1 is
        ## the spline only. the spline is constrained to its null-space
        phi.intern = list(
            old = c(),
            new = c()), 
        phi = list(
            old = c(),
            new = c()), 
        prec = list(
            old= c(),
            new = c()), 

        ## holds the complete covariance matrix
        S.full = list(
            old = c(),
            new = c()),

        lprior = list(
            old = c(),
            new = c())
        )
    )

for(k in 1:model$n.terms){
    xval = matrix(scale(model$X[[k]]), model$n, 1)
    model$S.fixed[[k]] = xval %*% t(xval) + model$delta.add.diag * diag(model$n)
    print(paste("compute S.fixed", k))

    values = get(paste(names(model$X)[k], ".values", sep=""))
    m = max(values)
    print(paste("random", k, "max.value=", m))
    A = matrix(0, model$n, m)
    for(i in 1:model$n) {
        A[i, model$X[[k]][i]] = 1
    }
    S.rw = INLA:::inla.ginv(INLA:::inla.rw(m, order = 2, scale.model=TRUE, sparse=FALSE))
    model$S.random[[k]] = A %*% S.rw %*% t(A) + model$delta.add.diag * diag(model$n)
    print(paste("compute S.random", k))
}

## build the three case 1 priors. they turn out to be almost identical...
for(k in 1:model$n.terms) {
    print(paste("compute priors.case1 ", k))
    model$priors.case1[[k]] = pcp.case1(
                          S.base = model$S.fixed[[k]],
                          S.alt = model$S.random[[k]],
                          U = model$param.case1$U,
                          alpha = model$param.case1$alpha,
                          use.phi.intern = TRUE)
    if (run.inla) {
        phi.intern = seq(-10, 10, len=1000)
        if (k==1) {
            inla.dev.new()
            plot(phi.intern, exp(model$priors.case1[[k]](phi.intern)), type="l", lty=k)
        } else {
            lines(phi.intern, exp(model$priors.case1[[k]](phi.intern)), lty=k, lwd=2*k)
        }
    }
}


## initialize 'old'
model$weights.intern$old = rnorm(model$n.terms -1L) ## yes, its one less...
model$weights$old = phi2w(model$weights.intern$old)  ## ...than here
model$phi.intern$old = rep(0, model$n.terms) ## one for each covariate
model$phi$old = inla.link.invlogit(model$phi.intern$old)
model$prec$old = 1 ## overall scaling

pc.joint.prior = function(weights.intern, phi.intern, log.prec, model)
{
    stopifnot(is.environment(model))
    
    model$weights.intern$old = weights.intern
    model$weights$old = phi2w(weights.intern) 
    model$phi.intern$old = phi.intern
    model$phi$old = inla.link.invlogit(phi.intern)
    model$prec$old = exp(log.prec)
    ##print(cbind(phi.intern = model$phi.intern$old,  phi = model$phi$old))

    S.all.covariates = list()
    for(k in 1:model$n.terms) {
        ##print("compute S.all.covariates", k)
        S.all.covariates[[k]] = ((1-model$phi$old[k]) * model$S.fixed[[k]] + model$phi$old[k] * model$S.random[[k]])
    }
    print(model$priors.case2$old)
    if (FALSE || length(model$priors.case2$old) == 0) {
        model$priors.case2$old = pcp.case2.H(S.all.covariates, eps=1e-3)
    } else {
        print("USE OLD case2 prior")
    }

    model$lprior$old = (inla.pc.dprec(prec = model$prec$old,
                                      u = model$param.prec$U,
                                      alpha = model$param.prec$alpha,
                                      log = TRUE) + log(model$prec$old)
                        +
                        inla.pc.multvar.sphere.general.d(
                            x = model$weights.intern$old,
                            lambda = model$param.case2$lambda,
                            log=TRUE,
                            H = model$priors.case2$old))
    for(k in 1:model$n.terms) {
        model$lprior$old = model$lprior$old +
            model$priors.case1[[k]](model$phi$old[k])
    }

    return (model$lprior$old)
}


inla2pc.mapping = function(x) 
{
    m = length(x) %/% 2
    lprec.beta = x[1:m]
    lprec.f = x[-(1:m)]
    stopifnot(length(lprec.beta) == length(lprec.f))
    
    ## mapping between INLA parameters and pc-prior parameters
    n = length(lprec.beta)
    prec.beta = exp(lprec.beta)
    prec.f = exp(lprec.f)
    phi = prec.beta/(prec.beta + prec.f)
    phi.intern = inla.link.logit(phi)

    w = numeric(n)
    w[1] = 1/sum(prec.beta[1]*(1-phi[1])/(prec.beta*(1-phi)))
    for(k in 2:n) {
        w[k] = w[1] * prec.beta[1]*(1-phi[1])/(prec.beta[k]*(1-phi[k]))
    }
    w.intern = w2phi(w)
    prec = prec.beta[1] * w[1] * (1-phi[1])
    x = list(weights.intern = w.intern, phi.intern = phi.intern, log.prec = log(prec))
    return(x)
}

log.jac = function(x)
{
    require(numDeriv)
    jac = jacobian(inla2pc.mapping, x)
    jac = unlist(jac)
    ## magic, this determinant is 1!!!!
    ##print(jac)
    ##print(det(jac))
    return (log(abs(det(jac))))
}
    

if (TRUE) {
    lprec.beta = rnorm(3)
    lprec.f = rnorm(3)
    x = unlist(inla2pc.mapping(c(lprec.beta, lprec.f)))
    w=phi2w(x[1:2])
    phi = inla.link.invlogit(x[3:5])
    prec = exp(x[6])
    print(cbind(log(prec/w/(1-phi)), lprec.beta))
    print(cbind(log(prec/w/(phi)), lprec.f))
    ##print(log.jac(x))
}

if (run.inla) {
    ## additive model in section 6.4
    initial = rnorm(6, sd = .1)
    r6.4 = inla(y ~ 1 +
        f(BP.idx, scale(BP), model="iid",
          hyper = list(prec = list(initial = initial[1]))) +
        f(CR.idx, scale(CR), model="iid", 
          hyper = list(prec = list(initial = initial[2]))) + 
        f(Age.idx, scale(Age), model="iid", 
          hyper = list(prec = list(initial = initial[3]))) +
        ## I have to reverse the order of CR and BP in order to get the hyperpar in inla in the
        ## correct order.
        f(CR,
          model="rw2",
          hyper = list(prec =
              list(prior = "pc.prec",
                   param = pc.prec.param,
                   initial = initial[5])), 
          values = CR.values, 
          constr = FALSE,
          diagonal = diag.value, 
          extraconstr = make.constr(CR.values)) + 
        f(BP,
          model="rw2",
          hyper = list(prec =
              list(prior = "pc.prec",
                   param = pc.prec.param,
                   initial = initial[4])), 
          values = BP.values, 
          constr = FALSE,
          diagonal = diag.value, 
          extraconstr = make.constr(BP.values)) + 
        f(Age,
          model="rw2",
          hyper = list(prec =
              list(prior = "pc.prec",
                   param = pc.prec.param,
                   initial = initial[6])), 
          values = Age.values, 
          constr = FALSE,
          diagonal = diag.value, 
          extraconstr = make.constr(Age.values)), 
        data = heart, 
        control.predictor = list(compute=TRUE),
        control.inla = list(tolerance = 1e-8, h = 0.001), 
        control.fixed = list(prec.intercept = 1/model$param.intercept$variance), 
        family = "binomial",
        inla.call = "./inla.special",
        verbose=TRUE)


    par(mfrow=c(3, 1))
    m = mean(heart$BP)
    s = sd(heart$BP)
    r = range(heart$sbp)
    plot(seq(r[1], r[2], len = length(BP.values)),
         r6.4$summary.random$BP.idx$mean[1] * (BP.values - m) /s +
         r6.4$summary.random$BP$mean, main = "BP")
    m = mean(heart$CR)
    s = sd(heart$CR)
    r = range(heart$ldl)
    plot(seq(r[1], r[2], len = length(CR.values)),
         r6.4$summary.random$CR.idx$mean[1] * (CR.values - m) /s +
         r6.4$summary.random$CR$mean, main = "CR")
    m = mean(heart$Age)
    s = sd(heart$Age)
    r = range(heart$age)     
    plot(seq(r[1], r[2], len = length(Age.values)),
             r6.4$summary.random$Age.idx$mean[1] * (Age.values - m) /s +
             r6.4$summary.random$Age$mean, main = "Age")

} else {

    fd1 = fifo("FIFO.FROM.INLA", "rb", blocking=TRUE)
    fd2 = fifo("FIFO.TO.INLA", "wb", blocking=TRUE)

    while (TRUE) {
        m = 2*model$n.terms
        print("Wait to read...")
        theta = readBin(fd1, what = double(), n = m)
        print(theta)
        my.theta = inla2pc.mapping(theta)
        ##print(my.theta)
        print(c(weights = phi2w(my.theta$weights.intern),
                phi = inla.link.invlogit(my.theta$phi.intern),
                log.prec = my.theta$log.prec))
        print("Evaluate joint pc-prior")
        lprior = pc.joint.prior(
            weights.intern = my.theta$weights.intern, 
            phi.intern = my.theta$phi.intern, 
            log.prec = my.theta$log.prec,
            model)
        print(paste("Done. lprior=", lprior))
        cat("\n\n")
        writeBin(as.double(lprior), fd2)
    }

}
