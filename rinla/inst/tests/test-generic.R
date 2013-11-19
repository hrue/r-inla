context("test 'latent generic'")

test_that("Case 1", {
    set.seed(1234)
    n = 100
    Q = INLA:::inla.rw1(n)
    i = 1:n
    y = (i-n/2)^2 + rnorm(n)
    
    r = inla(y ~ f(i, model="rw1", param = c(1, .001), constr=FALSE) - 1,
            data = data.frame(i,y))
    rr = inla(y ~ f(i, model="generic", Cmatrix=Q,
            rankdef=1, param=c(1, .001), constr=FALSE) - 1,
            data = data.frame(i,y))
    expect_true(all(abs(rr$summary.random[[1]] -
                        r$summary.random[[1]]) < 0.0001))
})



