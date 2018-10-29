########################################################################################
## Code for "rats" example, illustrating the use of the inla.cut()                    ##
## function for group-wise model criticism; see ?inla.cut for documentation.          ## 
##                                                                                    ##
## For details see Section 4.1 of the following paper:                                ##
## Egil Ferkingstad, Leonhard Held and Håvard Rue.                                    ##
## Fast and accurate Bayesian model criticism and conflict diagnostics using R-INLA.  ##
## Published in Stat, doi: 10.1002/sta4.163, 2017.                                    ##     
## Available as arXiv preprint arXiv:1708.03272: http://arxiv.org/abs/1708.03272      ##
########################################################################################


library(INLA)
## RATS DATA:
## Each row of ratsdata$y corresponds to a rat i, each column to a time point j.
## The value of ratsdata$y[i,j] is the weight of rat i at time j.
## See Ferkingstad et al. (2017) for further details about the data and model.
ratsdata <- structure(list(t = c(8, 15, 22, 29, 36), tbar = 22, N = 30, T = 5, 
    Omega = structure(c(200, 0, 0, 0.2), .Dim = c(2L, 2L)), mean = c(0, 
    0), prec = structure(c(1e-06, 0, 0, 1e-06), .Dim = c(2L, 
    2L)), y = structure(c(151, 145, 147, 155, 135, 159, 141, 
    159, 177, 134, 160, 143, 154, 171, 163, 160, 142, 156, 157, 
    152, 154, 139, 146, 157, 132, 160, 169, 157, 137, 153, 199, 
    199, 214, 200, 188, 210, 189, 201, 236, 182, 208, 188, 200, 
    221, 216, 207, 187, 203, 212, 203, 205, 190, 191, 211, 185, 
    207, 216, 205, 180, 200, 246, 249, 263, 237, 230, 252, 231, 
    248, 285, 220, 261, 220, 244, 270, 242, 248, 234, 243, 259, 
    246, 253, 225, 229, 250, 237, 257, 261, 248, 219, 244, 283, 
    293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 273, 289, 
    326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 286, 
    303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
    338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 
    321, 334, 302, 302, 323, 331, 345, 333, 316, 291, 324), .Dim = c(30L, 
    5L))), .Names = c("t", "tbar", "N", "T", "Omega", "mean", 
"prec", "y"))

## Specify the data in a format suitable for INLA:
inla_data <- data.frame(
  y = as.vector(ratsdata$y),
  tc = rep(ratsdata$t - mean(ratsdata$t), each=ratsdata$N), ## using centered t
  i1 = rep(1:ratsdata$N,ratsdata$T),                        ## rat index
  i2 = ratsdata$N + rep(1:ratsdata$N,ratsdata$T),           ## "copy-rat" index
  intercept = 1                                             ## explicit intercept 
)

## Formula defining model:
inla_formula <- 
  y ~ -1 + intercept + tc +
  f(i1, model="iid2d", n=2*30,
    hyper = list(theta1 = list(param = c(3, 200, .2, 0)))) +
  f(i2, tc, copy="i1")

## Initial INLA run for data and model:
r <- inla(formula = inla_formula,
               family = "gaussian",
               data = inla_data,
               control.predictor = list(compute=TRUE),
               control.family = 
                   list(hyper = list(theta = list(prior="loggamma", param=c(.001,.001)))),
               control.fixed = list(prec.intercept = 1E-6, prec = 1E-6),
               control.inla = list(h=0.03)
               )

## Run model criticism algorithm, giving a conflict p-value for each individual rat:
p <- inla.cut(r, split.by = "i1", debug=TRUE)

## Plot empirical distribution of p-values,
## with red line indicating false discovery rate threshold of 10%:
n.p <- length(p)
grid <- c(0:n.p)/n.p
y <- 0+0.1*grid
sel <- (sort(p) < y[-1])
par(pty="s")
plot(c(1:length(p))/n.p, sort(p), xlim=c(0,1), ylim=c(0,1), pch=19, cex=0.75,
     xlab="Rank of P-values divided by total number (30) of P-values", 
     ylab="Conflict P-values (ordered)", col=c(rep(2, sum(sel)), rep(1, n.p-sum(sel))))
abline(0,1, lty=2)
lines(grid, y, lty=2, col=2, type="l")
text(0.9, 0.15, "FDR = 10%", col=2)
abline(0,0.1, lty=1, col=2)

## Rat number 9 is seen to be divergent at FDR = 10%:
print(which.min(p))

## Conflict P-value for rat 9:
print(min(p))
