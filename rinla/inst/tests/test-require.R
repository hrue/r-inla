context("test 'require'")

test_that("INLA library load correctly", {
    expect_true(require(INLA))
})
