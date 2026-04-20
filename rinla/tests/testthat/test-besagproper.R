context("test 'latent besagproper'")

test_that("Case 1", {
    set.seed(12)
    ## pick a graph
    graph.file = system.file("demodata/germany.graph", package="INLA")
    g = inla.read.graph(graph.file)

    ## we will use replicated samples in our testing
    nrep = 5

    ## make life easy; use dense matrix algebra
    d = 1.0
    tau = 1.0
    Q = matrix(0, g$n, g$n)
    diag(Q) = tau * (d + g$nnbs)
    for(i in 1:g$n) {
        if (g$nnbs[i] > 0) {
            Q[i, g$nbs[[i]]] = -tau
            Q[g$nbs[[i]], i] = -tau
        }
    }
    R = chol(Q) ## 'chol' returns the upper triangular
    ## simulate data with replications
    y = c()
    for(i in 1:nrep) {
        y = c(y, backsolve(R, rnorm(g$n)))
    }
    i = rep(1:g$n, nrep)
    replicate = rep(1:nrep, each = g$n)
    formula = y ~ f(i, model="besagproper",  graph = graph.file,
            replicate=replicate,
            hyper = list(diag = list(param = c(1, 1)))) -1
    ## use 'exact' observations, so we fix the noise precisin to a high
    ## value
    r = inla(formula,
            data = data.frame(y, i, replicate),
            family = "gaussian",
            control.family = list(
                    hyper = list(
                            prec = list(
                                    initial = 10,
                                    fixed=TRUE))))
    expect_true(abs(r$summary.hyperpar[1, "mean"] < d) < 0.05)
})
