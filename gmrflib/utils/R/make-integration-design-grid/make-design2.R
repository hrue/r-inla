x <- c(0.0, 0.50, 1.00, 1.75, 2.50, 3.5)
x <- unique(sort(x))
m <- length(x)
w <- rep(1, m)

opt.fun <- function(par, args) {
    x <- args$x
    w <- args$w
    m <- args$m
    w[2:m] <- w[-1] * exp(par)
    ww <- c(rev(w[-1]), w)
    xx <- c(-rev(x[-1]), x)
    ww <- ww * dnorm(xx)
    ww <- ww/sum(ww)
    f2 <- sum(xx^2 * ww)
    f4 <- sum(xx^4 * ww)
    f6 <- sum(xx^6 * ww)
    f8 <- sum(xx^8 * ww) 
    f10 <- sum(xx^10 * ww) ## 945
    return ((f2-1)^2 + (f4-3)^2 + (f6-15)^2 + (f8-105)^2)
}

sol <- optim(rep(0, m-1), opt.fun,
             args = list(x = x, w = w, m = m),
             control = list(ndeps = rep(0.001, m-1)), 
             method = "BFGS")

w <- c(1, exp(sol$par))
ww <- c(rev(w[-1]), w)
xx <- c(-rev(x[-1]), x)

v <- ww * dnorm(xx)
v <- v / sum(v)
for(k in 1:8) {
    print(c(k, sum(xx^k * v)))
}

fnm <- "design-dim1.txt"
write(length(xx), file = fnm)
write(t(cbind(xx, ww)), file = fnm, append = TRUE, ncol = 2)

fnm <- "design-dim1-x.txt"
write(length(xx), file = fnm)
write(xx, file = fnm, append = TRUE, ncol = 1)

fnm <- "design-dim1-w.txt"
write(length(ww), file = fnm)
write(ww, file = fnm, append = TRUE, ncol = 1)


##############################################################
##############################################################
##############################################################

x <- c(0.0, 0.5, 1.25, 2.25)
x <- unique(sort(x))
m <- length(x)
w <- rep(1, m)

opt.fun <- function(par, args) {
    x <- args$x
    w <- args$w
    m <- args$m
    w[-1] <- w[-1] * exp(par)
    ww <- c(rev(w[-1]), w)
    xx <- c(-rev(x[-1]), x)
    ww <- ww * dnorm(xx)
    ww <- ww/sum(ww)
    f2 <- sum(xx^2 * ww)
    f4 <- sum(xx^4 * ww)
    ##f6 <- sum(xx^6 * ww)
    ##f8 <- sum(xx^8 * ww)
    ##f10 <- sum(xx^10 * ww)
    return ((f2-1)^2 + (f4-3)^2)
}

sol <- optim(rep(0, m-1), opt.fun,
             args = list(x = x, w = w, m = m),
             control = list(ndeps = rep(0.001, m-1)), 
             method = "BFGS")

w <- c(1, exp(sol$par))
ww <- c(rev(w[-1]), w)
xx <- c(-rev(x[-1]), x)

v <- ww * dnorm(xx)
v <- v/sum(v)
for(k in 1:8) {
    print(c(k, sum(xx^k * v)))
}

m <- length(xx)
vv <- matrix(NA, ncol = 3, nrow=m^2)
k <- 1
for(i in 1:m) {
    for(j in 1:m) {
        vv[k, ] <- c(xx[i], xx[j], ww[i] * ww[j])
        k <- k+1
    }
}

fnm <- "design-dim2.txt"
write(nrow(vv), file = fnm, append = FALSE)
write(t(vv), file = fnm, append = TRUE, ncol = 3)

fnm <- "design-dim2-x.txt"
write(nrow(vv), file = fnm, append = FALSE)
write(t(vv[, 1:2]), file = fnm, append = TRUE, ncol = 2)

fnm <- "design-dim2-w.txt"
write(nrow(vv), file = fnm, append = FALSE)
write(t(vv[, 3]), file = fnm, append = TRUE, ncol = 1)

if (FALSE) {
    ## check
    s0 <- 0
    s2 <- 0
    s3 <- 0
    s4 <- 0
    for(k in 1:m^2) {
        x1 <- vv[k, 1]
        x2 <- vv[k, 2]
        d <- dnorm(x1) * dnorm(x2) * vv[k, 3]
        s0 <- s2 + d * 1
        s2 <- s2 + d * x1^2
        s4 <- s4 + d * x1^4
    }
    print(c(variance = s2/s0, kurt = s4/s0))
}
