context("test 'expand-factor-strategy'")

test_that("Case 1", {
    r = inla(y ~ 1 + x + xx,
            data = data.frame(
                    y=1:4,
                    x=factor(c("a","b","c", NA)), 
                    xx=factor(c("A","B","C", NA))), 
            control.fixed = list(expand.factor.strategy="inla"))

    expect_true("xa" %in% colnames(r$model.matrix))
    expect_true("xb" %in% colnames(r$model.matrix))
    expect_true("xc" %in% colnames(r$model.matrix))
    expect_true("xxA" %in% colnames(r$model.matrix))
    expect_true("xxB" %in% colnames(r$model.matrix))
    expect_true("xxC" %in% colnames(r$model.matrix))
})

