## One-sample Kolmogorov-Smirnov test, plotting the normalised
## deviation between the empirical distribution function and the null
## hypothesis.
##
## Example:
## result = inla(..., control.compute=list(cpo=TRUE))
## inla.ks.plot(result$pit, punif)
inla.ks.plot = function (x, y, diff=TRUE, ...)
{
    if (any(is.na(x))) {
        x = x[!is.na(x)]
    }
    test = ks.test(x, y, ...)
    n = length(x)
    Fn = ((1:n)-0.5)/n
    F = y(sort(x))
    empirical.diff = (Fn-F)*sqrt(n)
    T = max(abs(empirical.diff))
    if (diff) {
        ylim = c(-1,1)*max(1,T)
        plot(F, empirical.diff, type='l',
             xlim=c(0,1),
             ylim=ylim,
             main=paste("K-S-test, p-value = ", test$p.value),
             ylab="(Fn-F) sqrt(n)",
             xlab="Quantile"
             )
        lines(c(0,1), c(0,0), type='l')
        lines(Fn, 2*sqrt(Fn*(1-Fn)), type='l')
        lines(Fn, -2*sqrt(Fn*(1-Fn)), type='l')
        lines(c(0,1,NA,0,1), c(T,T,NA,-T,-T), type='l')
    } else {
        plot(F, Fn, type='l',
             xlim=c(0,1),
             ylim=c(0,1),
             main=paste("K-S-test, p-value = ", test$p.value),
             ylab="Fn and F",
             xlab="Quantile"
             )
        lines(c(0,1), c(0,1), type='l')
        lines(Fn, Fn+2*sqrt(Fn*(1-Fn)/n), type='l')
        lines(Fn, Fn-2*sqrt(Fn*(1-Fn)/n), type='l')
        lines(c(0,1,NA,0,1), c(0,1,NA,0,1)+c(T,T,NA,-T,-T)/sqrt(n), type='l')
    }
    invisible(test)
}
