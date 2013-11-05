context("test 'latent besagproper2'")

test_that("Case 1", {
    set.seed(1234)
    ## pick a graph
    graph.file = system.file("demodata/germany.graph", package="INLA")
    g = inla.read.graph(graph.file)
    ## we will use replicated samples in our testing
    nrep = 5
    tau = 10.0
    lambda = 0.3
    R = -inla.graph2matrix(g)
    diag(R) = g$nnbs
    n = g$n
    Q = tau * ( (1-lambda) * diag(n) + lambda * R)
    y = c(inla.qsample(nrep, Q))
    i = rep(1:g$n, nrep)
    replicate = rep(1:nrep, each = g$n)
    formula = y ~ f(i, model="besagproper2",  graph = g,
            replicate=replicate) - 1
    r = inla(formula,
            data = data.frame(y, i, replicate),
            family = "gaussian",
            control.family = list(
                    hyper = list(
                            prec = list(
                                    initial = 10,
                                    fixed=TRUE))))
    expect_true(abs(r$summary.hyperpar[2, "mean"] - lambda) < 0.05)
})
