mi <- 0.33
mj <- 0.44
Sii <- 0.55
Sij <- 0.66
Sjj <- 0.88

S <- matrix(c(Sii, Sij, Sij, Sjj), 2, 2)
library(mvtnorm)

n <- 10^6
x <- rmvnorm(n, mean = c(mi, mj), sigma = S)

for(i in 1:10) {
    cat("moment(", i,  ") ",  mean(x[, 1]^i), "\n", sep = "")
}

ref <- readLines("ref")
k <- 1
for(i in 1:5) {
    for(j in 1:5) {
        cat("moment(", i, ", ", j, ") ",  mean(x[, 1]^i * x[, 2]^j), " ", as.numeric(ref[k]),
            "\n", sep = "")
        k <- k + 1
    }
}
