library(logspline)
burnin = 1:500
all.fnm = c(
        "result--lambda-0.1--exch-0.Rsave",
        "result--lambda-1--exch-0.Rsave",
        "result--lambda-0.1--exch-1.Rsave",
        "result--lambda-1--exch-1.Rsave")
new = TRUE
lty = 1
target = 3
for(fnm in all.fnm){
    load(fnm)
    ## beta2
    xx = result$beta.matrix[-burnin, target]
    if (new) {
        plot(logspline(xx),
             main = "",
             bty = "l", 
             ylim = c(0, 4), lty=lty, lwd=2,
             ylab = "Density",
             xlab = "Effect of smooking")
        new = FALSE
    } else {
        plot(logspline(xx),
             lwd = 2, lty=lty, add=TRUE)
    }
    print(paste("beta: column", target, fnm, mean(xx),  sd(xx),
                sum(xx>0)/length(xx)))
    lty = lty + 1
}
dev.print(postscript, file="mvprobit-beta2.ps")
system("which psfix && psfix mvprobit-beta2.ps")


load(all.fnm[3])
r = result$r.matrix
plot(logspline(r[, 1]), main = "",  bty = "l",
     lty = 1, lwd = 3, 
     xlim = c(0, 1),  xlab = "Correlation", ylab = "Density")
load(all.fnm[1])
r = result$r.matrix
lty = 1
for(k in 1:ncol(r)) {
    plot(logspline(r[, k]), add = TRUE,
         lwd = 1, lty = lty)
    lty = lty + 1
}
dev.print(postscript, file="mvprobit-r-lambda01.ps")
system("which psfix && psfix mvprobit-r-lambda01.ps")


load(all.fnm[4])
r = result$r.matrix
plot(logspline(r[, 1]), main = "",  bty = "l",
     lty = 1, lwd = 3, 
     xlim = c(0, 1),  xlab = "Correlation", ylab = "Density")
load(all.fnm[2])
r = result$r.matrix
lty = 1
for(k in 1:ncol(r)) {
    plot(logspline(r[, k]), add = TRUE,
         lwd = 1, lty = lty)
    lty = lty + 1
}
dev.print(postscript, file="mvprobit-r-lambda1.ps")
system("which psfix && psfix mvprobit-r-lambda1.ps")



