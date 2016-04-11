make.new.figures = TRUE

rarp = function(p, lambda, ...)
{
    pac = inla.pc.ar.rpacf(n=1, p, lambda)
    phi = inla.ar.pacf2phi(pac)
    return (phi)
}

lambda = 0.5
nsim = 100000
p=4
require(parallel)
X = mclapply(1:nsim, rarp, p=p, lambda=lambda, mc.cores = detectCores())
X = matrix(unlist(X),ncol=p, byrow=TRUE)
XX = apply(X, 1, function(x) inla.ar.phi2acf(x)[-1])
XX = matrix(unlist(XX), ncol = p,  byrow = TRUE)
colnames(X) = c("psi[1]", "psi[2]", "psi[3]", "psi[4]")
pairs(X[1:(nsim/10), ], cex=0.1,
      labels = c(expression(psi[1]), expression(psi[2]),
              expression(psi[3]), expression(psi[4])))
              
if (make.new.figures) {
    dev.print(postscript, file="ar4-joint.ps")
}

if (TRUE) {
    ##
    inla.dev.new()
    par(usr=c(-2000, 2000, -2000, 2000))
    par(mfrow = c(4, 4))
    nsim.plot = nsim %/% 10L
    for(i in 1:p) {
        for(j in 1:p) {
            print(c(i, j))
            if (i == j) {
                xx = c(X[, i], abs(X[, i]))
                xx = c(xx, -xx)
                rr = max(xx)
                hist(xx, prob = TRUE,
                     xlim = c(-rr, rr),
                     breaks = seq(-rr, rr, len=50*p),
                     xlab = expression(psi [1]),
                     ylab = "Density",
                     main = "")
            } else if (i < j) {
                plot(X[1:nsim.plot, j], X[1:nsim.plot, i],
                     cex=0.1,
                     bty="l", 
                     xlab = paste("phi[", j, "]", sep=""),
                     ylab = paste("phi[", i, "]", sep=""), 
                     xlim = range(X[, j]), ylim=range(X[, i]))
            } else if (j < i) {
                plot(XX[1:nsim.plot, j], XX[1:nsim.plot, i],
                     cex=0.1,
                     bty="l", 
                     xlab = paste("acf[", j, "]", sep=""),
                     ylab = paste("acf[", i, "]", sep=""), 
                     xlim = c(-1, 1), ylim = c(-1, 1))
            } else {
                stop("SHN")
            }
        }
    }
    if (make.new.figures) {
        dev.print(postscript, file="ar4-all.ps")
    }
}


## display the marginal for psi[1]

xx = X[, 1]
xx = abs(xx)
xx = c(xx, -xx)
inla.dev.new()
hist(xx, prob = TRUE,  xlim = c(-p, p),  breaks = seq(-(p+1), p+1, len=50*p),
     xlab = expression(psi [1]),
     ylab = "Density",
     main = "")
if (make.new.figures) {
    dev.print(postscript, file="ar4-phi1.ps")
}

if (make.new.figures) {
    system("which psfix && psfix ar4-*.ps")
}


