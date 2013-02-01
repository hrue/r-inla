n = 300
beta = 2
x = rnorm(n)
prec.x = 10
prec.y = 1000
s = runif(n)
x.tilde = x + rnorm(n, sd = 1/sqrt(s*prec.x))
y = 1 + beta * x.tilde + rnorm(n, sd = 1/sqrt(prec.y))

r = inla(y ~ f(x, model="meb", scale = s), 
        family = "gaussian", 
        data = data.frame(y, x, s))
