n <- 200
loc <- seq(0, n-1, length.out=n)
h.size <- diff(range(loc)) / (n-1)
f.true <- (sin(2*pi*(loc/n)^2))^3
s = 0.2
y <- f.true + rnorm(n, sd = s)

plot(loc, y, pch = 19)
lines(loc, f.true, type='l', lwd = 3, col='blue')

r <- inla(y ~ -1 + f(loc, model = "prw2", values = loc,
                     hyper = list(
                         range = list(
                             param = c(50, 0.5, h.size, 0)))), 
          family = "normal",
          control.family = list(
              hyper = list(prec = list(initial = log(1/s^2),
                                       fixed = TRUE))), 
          data = data.frame(y, loc))
lines(loc, r$summary.linear.predictor$mean, lwd = 3, col = "red")

## we can now change the resolution and get the 'same' results
nn <- 2*n
lloc <- seq(0, n-1, length.out = nn)
yy <- numeric(nn)
for(i in 1:n) {
    yy[1 + (i-1) * 2] <- y[i]
    yy[2 + (i-1) * 2] <- NA
}
hh.size <- diff(range(lloc)) / (nn-1)

rr <- inla(yy ~ -1 + f(lloc, model = "prw2", values = lloc,
                     hyper = list(
                         ## note that '50' is the same as its in the real scale
                         range = list(param = c(50, 0.5, hh.size, 0)))), 
          family = "normal",
          control.family = list(hyper = list(prec = list(initial = log(1/s^2),
                                                         fixed = TRUE))), 
          data = data.frame(yy, lloc))
lines(lloc, rr$summary.linear.predictor$mean, lwd = 3, col = "red")

## once more
nnn <- 4*n
llloc <- seq(0, n-1, length.out = nnn)
yyy <- numeric(nnn)
for(i in 1:n) {
    yyy[1 + (i-1) * 4] <- y[i]
    yyy[2:4 + (i-1) * 4] <- NA
}
hhh.size <- diff(range(llloc)) / (nnn-1)

rrr <- inla(yyy ~ -1 + f(llloc, model = "prw2", values = llloc,
                     hyper = list(
                         ## again, '50' is not changed
                         range = list(param = c(50, 0.5, hhh.size, 0)))), 
          family = "normal",
          control.family = list(hyper = list(prec = list(initial = log(1/s^2),
                                                         fixed = TRUE))), 
          data = data.frame(yy, lloc))
lines(llloc, rrr$summary.linear.predictor$mean, lwd = 3, col = "red")

inla.dev.new()
res <- cbind(r$summary.linear.predictor$mean,
             rr$summary.linear.predictor$mean[seq(1, nn, by = 2)], 
             rrr$summary.linear.predictor$mean[seq(1, nnn, by = 4)])
pairs(res)

inla.dev.new()
plot(r$marginals.hyperpar[[2]],  pch = 19)
lines(rr$marginals.hyperpar[[2]], lwd = 3, col = "blue")
lines(rrr$marginals.hyperpar[[2]], lwd = 3, col = "red")
