ll <- function(x, a) {
    phi <- 1/(1+exp(-x))
    return (log (a * phi + (1-phi)))
}

ll.new <- function(x, a) {
    phi <- 1/(1+exp(-x))
    x0 <- (-0.5 * log(a)) * 0.9
    x1 <- (-0.5 * log(a)) * 0.98
    idx <- which(x >= x0)
    if (length(idx) > 0) {
        xx <- x[idx]
        x0 <- x0 + (x1-x0) * (1 - exp(-(xx-x0)))
        ex0 <- exp(-x0)
        yy <- log((a+ex0) / (1+ex0)) +
            ((ex0*(a-1)) / ((a + ex0)*(1+ex0)))*(xx-x0) +
            0.5 * (((-a + (ex0)^2)*(a-1)*ex0) / ( (a+ex0)^2 * (1+ex0)^2)) * (xx-x0)^2
    }
    val <- log (a * phi + (1-phi))
    if (length(idx) >0) {
        val[idx] <- yy
    }
    return (val)
}

for (a in 1/10^(4)) {

    x <- seq(-10, 6, len = 10000)
    y <- ll(x, a)
    fun <- splinefun(x, y)
    ynew <- ll.new(x, a)
    fun.new <- splinefun(x, ynew)

    dev.new()
    par(mfrow = c(1, 2))
    plot(x,y, type = "l")
    lines(x,ynew, type = "l")
    plot(x,fun(x, deriv = 2), type = "l")
    lines(x,fun.new(x, deriv = 2), type = "l")
    abline(h = 0)
    title(paste0("a=", a))
}
