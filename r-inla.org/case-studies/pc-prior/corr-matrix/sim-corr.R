rcorr = function(p, lambda)
{
    theta = inla.pc.cormat.rtheta(n=1, p, lambda)
    R = inla.pc.cormat.theta2R(theta)
    R = inla.pc.cormat.permute(R)
    return (R)
}

make.new.figures = TRUE
first.time = TRUE
lty=1
nsim = 10000
p=3
m = p*(p-1)/2
X = matrix(NA, nsim, m)
for (lambda in c(10, 5, 2)) {
    for(k in 1:nsim) {
        r = rcorr(p=p, lambda = lambda)
        X[k, ] = r[upper.tri(r)]
    }
    if (p==3) colnames(X) = c("r12", "r13",  "r23")

    require(logspline)
    X = abs(X)
    d = logspline(c(X[, 1], -X[, 1]), lbound=-1, ubound=1)
    eps = 1e-5
    xx = seq(0, 1-eps, len=1000)
    xx = unique(sort(c(-xx, xx)))
    if (first.time) {
        first.time = FALSE
        plot(xx, dlogspline(xx, d), 
             lwd=2,
             lty=lty,
             type="l",
             bty = "l", 
             xlim = c(-1, 1),
             ylim = c(0, 12), 
             ylab = "Density",
             xlab = "Correlation")
    } else {
        lines(xx, dlogspline(xx, d), lwd=2, lty=lty)
    }
    lty = lty + 1
}

if (make.new.figures) {
    dev.print(postscript, file="corrmat3.ps")
    system("which psfix && psfix corrmat*.ps")
}

