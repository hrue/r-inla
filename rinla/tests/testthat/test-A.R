context("test 'the A-matrix'")

test_that("Case 1", {
    set.seed(123)
    eps = 0.025
    n = 300
    z = rnorm(n)
    zz = 1:n
    y = z + zz + rnorm(n, sd = 0.1)
    r = sample(1:n, n)
    Y = y[r]
    A=sparseMatrix(i=1:n, j=r, x=rep(1,n))
    AA = as.matrix(A)
    formula = Y ~ -1 + z + zz
    r = inla(formula, data = data.frame(Y,z,zz),
            control.predictor = list(compute=T, A = A))
    rr = inla(formula, data = data.frame(Y,z,zz),
            control.predictor = list(compute=T, A = AA))
    
    expect_true(all(abs(r$summary.fixed[, "mean"] - 1) < eps))
    expect_true(all(abs(rr$summary.fixed[, "mean"] - 1) < eps))
    expect_true(all(abs(r$summary.fixed- rr$summary.fixed) < eps))
})

test_that("Case 2", {
    set.seed(1234)
    eps = 0.05
    n = 100
    z = rnorm(n)
    zz = 1:n
    lin.pred = z + zz
    y = numeric(n)
    lin.pred.new = numeric(n)
    A=list(i=c(), j=c(), values=c())
    for(i in 1:n) {
        lin.pred.new[i] = mean(lin.pred[1:i])
        y[i] = lin.pred.new[i] + rnorm(1, sd=0.1)
        A$i = c(A$i, rep(i, i))
        A$j = c(A$j, 1:i)
        A$values = c(A$values, rep(1/i, i))
    }
    Am = sparseMatrix(i = A$i, j = A$j, x = A$values)
    formula = y ~ -1 + z + zz
    r = inla(formula, data = data.frame(y,z,zz),
            control.predictor = list(compute=T, A = Am))
    expect_true(all(abs(r$summary.linear.predictor$mean[1:n] - lin.pred.new) < eps))
    expect_true(all(abs(r$summary.linear.predictor$mean[1:n + n] -  lin.pred) < eps))
})


test_that("Case 3", {
    set.seed(12345)
    s = 0.00001
    n = 100
    fac = 2
    m = n*fac
    z = runif(n)
    zz = rep(z, fac)
    
    Y = rep(NA, m)
    A = matrix(0,m,m)
    for(i in 1:m) {
        ## just make some random selection
        A[i, sample(1:m, n %/% 10)] = 1
        Y[i] =  sum(A[i,]*(1+zz)) + rnorm(1, sd=s)
    }
    rr = inla(Y ~ -1 + zz, control.predictor = list(compute=TRUE,A=A),
            data = data.frame(z,Y),
            control.family = list(initial = log(1/s^2), fixed=T))
    expect_true(all(abs(rr$summary.linear.predictor$mean[1:n] - Y[1:n]) < sqrt(s)))
})

test_that("Case 4", {
    set.seed(123)
    n = 30
    m = 10
    off = 10+runif(m)
    z = runif(n)
    eta = 0 + z
    A = matrix(runif(n*m),m,n)
    s = 0.01
    Eta = A %*% eta + off
    Y = Eta +  rnorm(m, sd=s)
    
    r = inla(Y ~ -1+z, control.predictor = list(A=A,compute=TRUE,
                                            precision = 1e6),
            data = list(Y=Y, z=z, off = off),
            offset = off, 
            control.compute=list(cpo=T),
            control.family = list(initial = log(1/s^2), fixed=TRUE))
    expect_true(all(abs(Eta - r$summary.linear.predictor[1:m,"mean"]) < sqrt(s)))
    expect_true(all(abs(eta - r$summary.linear.predictor[m+1:n,"mean"]) < sqrt(s)))
})
