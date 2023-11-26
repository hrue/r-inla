n <- 1000
x <- rnorm(n)
xx <- rnorm(n)
off <- runif(n)
z <- rnorm(n) 
zz <- rnorm(n)

mean <- off + 0.1 + 1.1 * z + 2.2 * zz
prec <- exp(1 + 0.55 * x + 1.1 * xx)

y <- mean + (1/sqrt(prec)) * rnorm(n)
Y <- inla.mdata(y, off, 1, z, zz)
r <- inla(Y ~ 1 + x + xx,
          data = list(Y = Y, off = off, x = x, xx = xx, z = z, zz = zz),
          family = "ggaussianS")
summary(r)
