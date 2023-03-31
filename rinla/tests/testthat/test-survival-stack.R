context("test 'survival through stack'")

test_hat("Case 1: weibullsurv", {

    data(Leuk)
    Leuk$day <- Leuk$time / 365

    ff1 <- inla.surv(day, cens) ~
        age + sex + wbc + tpi
    ff1
    fit1 <- inla(
        formula = ff1,
        family = "weibullsurv",
        data = Leuk)
    fit1$mode$theta

    dstack <- inla.stack(
        data = list(
            y = Leuk$day,
            ce = Leuk$cens),
        effects = list(
            data.frame(
                a0 = 1,
                Leuk[c("age", "sex", "wbc", "tpi")])),
        A = list(1))

    ff2 <- inla.surv(y, ce) ~
        0 + a0 + age + sex + wbc + tpi
    ff2
    fit2 <- inla(
        formula = ff2,
        family = "weibullsurv",
        data = inla.stack.data(dstack))

    expect_true(abs(fit1$mode$theta -
                    fit2$mode$theta) < 0.1)

    sv <- inla.surv(Leuk$day, Leuk$cens)
    class(sv)
    str(sv)

    ## hack 1: to work for the names
    attr(sv, 'class') <- c("inla.surv", "data.frame")
    is.data.frame(sv)

    ## This is weird as the vectors in the output of `inla.surv()`
    ##  does not have all the same length.
    ## This may be solved if the inla.stack() methods consider
    ##  first `list` then `data.frame`:
    ##  then we only need as.list.inla.surv() 
    ##  and will not need the next hack.

    ## hack 2: to properly extract the data
    ## (weird to have to do it if it "is" `data.frame` already)
    as.data.frame.inla.surv <-
        function(x, ...) {
            class(x) <- "list"
            n <- length(x$time)
            nn <- lapply(x, length)
            return(as.data.frame(
                x[nn == n], ### is this always OK???
                ...))
        }

    dstackS <- inla.stack(
        data = list(sv[1:5]), 
        effects = list(
            data.frame(
                a0 = 1,
                Leuk[c("age", "sex", "wbc", "tpi")])),
        A = list(1))

    ## TO DO: make inla.stack() to collect the "names.ori"
    str(inla.stack.data(dstackS))
        
    ## problem: how do user use this to fit the model ???
    ## for now the problem is that we have to use the
    ## outcome names from inla.surv()

    ## the outcome names now are determined by inla.surv()
    ##   solution (proposed with Denis, useful for his code as well):
    ## add this in the FIRST line of inla.surv() code:
###       names.ori <- as.list(match.call())[-1]
    ## and add this BEFORE the LAST line of inla.surv() code:
###       attr(ret, "names.ori") <- names.ori
    
})


test_hat("Case 2: coxph", {

    data(Leuk)
    Leuk$day <- Leuk$time / 365
    ff1 <- inla.surv(day, cens) ~
        age + sex + wbc + tpi
    ff1
    fit1 <- inla(
        formula = ff1,
        family = "coxph",
        data = Leuk)
    fit1$mode$theta

    dstack <- inla.stack(
        data = list(
            y = Leuk$day,
            ce = Leuk$cens),
        effects = list(
            data.frame(
                a0 = 1,
                Leuk[c("age", "sex", "wbc", "tpi")])),
        A = list(1))

    ff2 <- inla.surv(y, ce) ~
        0 + a0 + age + sex + wbc + tpi
    ff2
    fit2 <- inla(
        formula = ff2,
        family = "coxph",
        data = inla.stack.data(dstack))

    expect_true(abs(fit1$mode$theta -
                    fit2$mode$theta) < 0.1)

})

