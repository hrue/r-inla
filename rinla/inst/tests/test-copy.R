context("test 'copy'")

test_that("Case 1", {
    set.seed(123)
    n = 300
    rho = 0.9
    sd = 0.01
    y = arima.sim(n=n, model = list(ar = rho))*sqrt(1-rho^2) + rnorm(n, sd = sd)
    y[1:10]=NA  ## ``burn in''
    i = 1:n
    j = i
    ii = c(NA, 1:(n-1))
    ww = rep(-1,n)
    na = rep(NA,n)
    
    i2 = c(i,i)
    ii2 = c(na, ii)
    j2 = c(na, j)
    ww2 = c(na, ww)
    
    Y = matrix(NA, 2*n, 2)
    Y[1:n, 1] = y
    Y[1:n + n, 2] = 0
    
    formula = Y ~ -1 + f(i2, initial=-10, fixed=T) +
        f(ii2, ww2, copy = "i2", range = c(-1,1),
          fixed = FALSE, initial = 3, precision = 1e10,
          param = c(0,0.2), prior = "normal") +
              f(j2, initial=0, fixed=FALSE)
    
    r = inla(formula, data = list(Y=Y, i2=i2,ii2=ii2,ww2=ww2,j2=j2),
            family = c("gaussian","gaussian"),
            control.family = list(
                    list(initial = log(1/sd^2), fixed = TRUE),
                    list(initial = 13, fixed = TRUE)),
            control.predictor = list(compute=TRUE))
    
    fformula = y ~ -1 + f(i, model = "ar1")
    rr = inla(fformula, data = data.frame(i,y),
            family = "gaussian",
            control.predictor = list(compute=TRUE),
            control.family = list(initial = log(1/sd^2), fixed = TRUE))
    h = inla.hyperpar(r,  diff.logdens = 10)
    hh = inla.hyperpar(rr,  diff.logdens=10)
    
    expect_true(abs(h$summary.hyperpar[2,"mean"] -
                    hh$summary.hyperpar[2,"mean"]) < 0.02)
})


test_that("Case 2", {
    set.seed(123)
    n = 500
    N = 4*n
    time = 1:n
    x = sin(time/n*4*pi)
    age = rep(1:4, each=n)
    
    y = numeric(N)
    y[ which(age==1) ] = x
    y[ which(age==2) ] = 2*x
    y[ which(age==3) ] = 3*x
    y[ which(age==4) ] = 4*x
    
    s = 0.01
    y = y + rnorm(N, sd=s)
    
    formula = y ~ -1 + f(i, model="rw2") +
        f(j, copy="i",  fixed=FALSE) + 
            f(k, copy="i",  fixed=FALSE) + 
                f(l, copy="i",  fixed=FALSE)
    
    i = rep(NA, N)
    i[which(age == 1)] = time
    j = rep(NA, N)
    j[which(age == 2)] = time
    k = rep(NA, N)
    k[which(age == 3)] = time
    l = rep(NA, N)
    l[which(age == 4)] = time
    r = inla(formula,  data = data.frame(y, i, j, k, l),  family = "gaussian",
            control.family = list(hyper =
                    list(prec = list(initial=log(1/s^2), fixed=TRUE))))
    expect_true(abs(r$summary.hyperpar[2, "mean"] - 2) < 0.05)
    expect_true(abs(r$summary.hyperpar[3, "mean"] - 3) < 0.05)
    expect_true(abs(r$summary.hyperpar[4, "mean"] - 4) < 0.05)
}) 
