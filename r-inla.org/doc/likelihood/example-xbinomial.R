n <- 10000
x <- rnorm(n, sd = 1)
q <- runif(n)
eta <- 0.88 + 0.77*x
p <- q * 1.0/(1+exp(-eta))
ntrials <- sample(1:25,  size=n, replace=TRUE)
y <- rbinom(n = n, size=ntrials, prob = p)
r <- inla(y ~ 1 + x,
         family = "xbinomial",
         Ntrials = ntrials,
         scale = q, 
         data = data.frame(y, x, q, ntrials))
summary(r)
