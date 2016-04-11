## plot all results...

burnin = 1:2000
all.fnm = dir(pattern="result--.*.Rsave")
load(all.fnm[1])
nbeta = ncol(result$beta.matrix)
nr = ncol(result$r.matrix)
for(i in 1:nbeta) {
    dev.new()
    new = TRUE
    lty = 1
    for(fnm in all.fnm) {
        load(fnm)
        xx = result$beta.matrix[-burnin, i]
        if (new) {
            plot(logspline(xx),
                 main = paste("beta", i, sep=""),
                 ylim = c(0, 8), lty=lty, lwd=2,
                 ylab = "Density")
            new = FALSE
        } else {
            plot(logspline(xx),
                 lwd = 2, lty=lty, add=TRUE)
        }
        print(paste("beta", i, fnm, mean(xx),  sd(xx)))
        lty = lty + 1
    }
}
stopifnot(nr == 6)
for(i in 1:nr) {
    dev.new()
    new = TRUE
    lty = 1
    for(fnm in all.fnm) {
        load(fnm)
        if (new) {
            plot(logspline(result$r.matrix[-burnin, i]),
                 main = paste("r", i, sep=""),
                 ylim = c(0, 8),
                 lty=lty, lwd=2,
                 ylab = "Density")
            new = FALSE
        } else {
            ii = min(i, ncol(result$r.matrix))
            plot(logspline(result$r.matrix[-burnin, ii]),
                 lwd = 2, lty=lty, add=TRUE)
        }
        lty = lty + 1
    }
}



