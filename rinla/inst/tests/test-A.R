context("testing 'the A-matrix'")

test_that("Case 1", {
    set.seed(123)
    eps = 0.05
    n = 100
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
    expect_true(all(abs(r$summary.fixed- rr$summary.fixed) < 0.05))
})
