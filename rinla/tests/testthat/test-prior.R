context("test 'prior spesification (table, expression)'")

test_that("Case 1", {
    set.seed(123)
    y = rnorm(100)
    a = 1
    b = 0.1
    hyper = list(prec = list(prior = "loggamma", param = c(a, b)))
    r = inla(y ~ 1,  data = data.frame(y),
            control.family = list(hyper = hyper))
    loggamma = "expression:
            a = 1;
            b = 0.1;
            precision = exp(log_precision);
            logdens = log(b^a) - lgamma(a)
                      + (a-1)*log_precision - b*precision;
            log_jacobian = log_precision;
            return(logdens + log_jacobian);"
    hyper.new = list(prec = list(prior = loggamma))
    r.new = inla(y ~ 1,  data = data.frame(y),
            control.family = list(hyper = hyper.new), verbose=FALSE)
    expect_true(abs(r$summary.hyperpar[1,"mean"] -
                    r.new$summary.hyperpar[1,"mean"]) < 0.001)
})

test_that("Case 2", {
    prior.function = function(log_precision) {
        a = 1;
        b = 0.1;
        precision = exp(log_precision);
        logdens = log(b^a) - lgamma(a) + (a-1)*log_precision - b*precision;
        log_jacobian = log_precision;
        return(logdens + log_jacobian)
    }
    prior.expression = "expression:
            a = 1;
            b = 0.1;
            precision = exp(log_precision);
            logdens = log(b^a) - lgamma(a)
                      + (a-1)*log_precision - b*precision;
            log_jacobian = log_precision;
            return(logdens + log_jacobian);"
          
    lprec = seq(-10, 10, len=1000)
    prior.table = paste(c("table:", cbind(lprec, prior.function(lprec))),
            sep = "", collapse = " ")
          
    set.seed(123)
    n = 100
    y = rnorm(n)
    r = inla(y~1,
            data = data.frame(y),
            control.family = list(
                    hyper = list(
                            prec = list(
                                    prior = "loggamma",
                                    param = c(1, 0.1)))))
    rr = inla(y~1,
            data = data.frame(y),
            control.family = list(
                    hyper = list(
                            prec = list(
                                    prior = prior.expression))))
    rrr = inla(y~1,
            data = data.frame(y),
            control.family = list(
                    hyper = list(
                            prec = list(
                                    prior = prior.table))))
    expect_true(abs(r$mlik[1] -  rr$mlik[1]) < 0.001)
    expect_true(abs(r$mlik[1] -  rrr$mlik[1]) < 0.001)
    expect_true(abs(rr$mlik[1] -  rrr$mlik[1]) < 0.001)
})
