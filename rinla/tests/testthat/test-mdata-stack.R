context("test 'mdata through stack'")

test_that("Case 1: 0poisson", {

    sim.poisson <- function(prob, m) {
        stopifnot(length(prob) == length(m) && length(prob) > 0)
        n <- length(m)
        y <- numeric(n)
        event <- (runif(n) < prob)
        idx.zero <- which(event)
        idx.non.zero <- which(!event)
        y[idx.zero] <- 0
        y[idx.non.zero] <- rpois(length(idx.non.zero), lambda = m[idx.non.zero])
        return (y)
    }

    n <- 300
    xx <- cbind(runif(n), runif(n))

    beta1 <- c(-1, 1, 0)
    prob <- plogis(cbind(1, xx) %*% beta1)

    beta2 <- c(1, 0, -1)
    lambda <- exp(cbind(1, xx) %*% beta2)
    
    E <- rgamma(n, 20, 2)

    repeat {
        if(any((y <- sim.poisson(prob, E*lambda))>0))
            break
    }

    dataf <- data.frame(
        y = y, E = E,
        x1 = xx[, 1], x2 = xx[, 2],
        x1y = xx[, 1], x2y = xx[, 2])

    ff1 <- inla.mdata(cbind(y, E),
                      cbind(1, x1, x2)) ~ 1 + x1y + x2y
    fit1 <- inla(
        formula = ff1,
        family = "0poisson",
        data = dataf)
    
    cbind(true = beta1, fit1$summary.hy[, 1:2])
    cbind(true = beta2, fit1$summary.fix[, 1:2])

### allow mdata going to inla.stack()
    
    mdtest <- inla.mdata(
        cbind(dataf$y, dataf$E),
        cbind(1, dataf$x1, dataf$x2))
    # str(mdtest)
    
    if(!any(attr(mdtest, "class") == "list")) {
        attr(mdtest, "class") <- c("inla.mdata", "list")
    }
        
    
    datas <- inla.stack(
        data = list(dummy = rep(NA, NROW(mdtest[[1]]))),
        response = list(mdtest), 
        effects = list(
            data.frame(Intercept = 1,
                       a0y = 1,
                       x1y = dataf$x1y,
                       x2y = dataf$x2y)),
        A = list(1))

    ff2 <- RESPONSE ~ 0 + Intercept + x1y + x2y
    fit2 <- inla(
        formula = ff2,
        data = inla.stack.data(datas, .response.name = "RESPONSE"),
        control.predictor = list(A = inla.stack.A(datas)),
        family = "0poisson")

    expect_equal(
        fit2$summary.hyperpar[, 1:2],
        fit1$summary.hyperpar[, 1:2],
        tolerance = 1e-5
    )
    rownames(fit2$summary.fixed)[1] <- "(Intercept)"
    expect_equal(
        fit2$summary.fixed[, 1:2],
        fit1$summary.fixed[, 1:2],
        tolerance = 1e-5
    )
})
