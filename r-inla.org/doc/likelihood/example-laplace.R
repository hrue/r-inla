

library(INLA)
inla.update(testing=T)

rald <- function(n, tau=0.5, mu=0, delta=1){
	xi <- (1-2*tau)/(tau - tau^2)		
	sigma <- sqrt(2/(tau - tau^2))
	w <- rexp(n, delta)
	mu + xi*w + sigma*sqrt(w/delta)*rnorm(n)
}

n <- 1000
mu <- rnorm(1)
tau <- 0.8
delta <- 0.5
y <- rald(n, tau, mu, delta)

formula <- y ~ mu - 1
r <- inla(formula, family = "laplace", data = data.frame(mu,y), 
			control.family = list(alpha=tau, gamma=2, epsilon=0.01))
summary(r)
