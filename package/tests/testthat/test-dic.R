context("test 'DIC'")

test_that("Case 1", {
    set.seed(123)
    n = 20
    y = rnorm(n)
    yy = rnorm(n)
    k=3
    y[1:k] = NA
    yy[n:n] = NA
    
    Y = matrix(NA, 2*n, 2)
    Y[1:n, 1] = y
    Y[n + 1:n, 2] = yy
    
    r = inla(Y ~1,  data = list(Y=Y),  family = rep("gaussian", 2),
            control.compute = list(dic = TRUE))

    expect_true(all(r$dic$family[(k+1):n] == 1))
    expect_true(all(is.na(r$dic$family[1:k])))
    expect_true(all(r$dic$family[n + 1:(n-1)] == 2))
    expect_true(all(is.na(r$dic$family[2*n])))

    expect_true(all(is.na(r$dic$local.dic[1:k])))
    expect_true(all(is.na(r$dic$local.p.eff[2*n])))
})



