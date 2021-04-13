library(sn)

power.13 <- function(x) sign(x) * abs(x)^(1/3)

find.sn.mode <- function(skew, mode.initial = NULL, plot = FALSE) {
    dp <- unlist(INLA:::inla.sn.reparam(moments = c(0, 1, skew)))
    if (is.null(mode.initial)) {
        mode.initial <- qsn(0.5, dp = dp)
    }
    res <- optim(mode.initial, dsn, log = TRUE, dp = dp,
                 method = "BFGS",
                 control = list(fnscale = -1, ndeps = 1e-5))

    if (plot) {
        xx <- seq(-3, 3, by = 0.001)
        my.plot(xx, dsn(xx, dp = dp), type = "l")
        abline(v = res$par, lwd = 3)
    }

    return(res$par)
}

find.d3.mode <- function(skew, mode) {
    dp <- unlist(INLA:::inla.sn.reparam(moments = c(0, 1, skew)))
    wfff = c(-7.0 / 240.0, 3.0 / 10.0, -169.0 / 120.0, 61.0 / 30.0, 0.0, -61.0 / 30.0, 169.0 /120.0, -3.0 / 10.0, 7.0 / 240.0 )
    dm<- 0.001
    ns <- length(wfff)
    ns2 <- ns %/% 2L + 1
    mm <- ns-ns2
    ms <- mode + seq(-mm, mm) * dm
    d3 <- sum(wfff * dsn(ms, dp = dp, log = TRUE) / dm^3)
    return(d3)
}


ss <- seq(-0.988, 0.988, by = 0.005)
mm <- numeric(length(ss))
d3 <- numeric(length(ss))

for(i in seq_along(ss)) {
    mm[i] <- find.sn.mode(skew = ss[i], mode.initial = if (i > 1) mm[i-1] else NULL)
    d3[i] <- find.d3.mode(skew = ss[i], mode = mm[i])
}

M <- matrix(ss, ncol = 1)
write(t(M), file = "skew-mode-matrix-s.txt", ncol = 1)
M <- matrix(power.13(ss), ncol = 1)
write(t(M), file = "skew-mode-matrix-s3.txt", ncol = 1)
M <- matrix(mm, ncol = 1)
write(t(M), file = "skew-mode-matrix-m.txt", ncol = 1)
M <- matrix(power.13(d3), ncol = 1)
write(t(M), file = "skew-mode-matrix-d33.txt", ncol = 1)
