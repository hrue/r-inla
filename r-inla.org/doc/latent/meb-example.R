n = 100
beta = 2
w = rnorm(n)
prec.u = 100
prec.y = 1000
s = runif(n, min = 0.5, max = 1/0.5)
x = w + rnorm(n, sd = 1/sqrt(s*prec.u))
y = 1 + beta * x + rnorm(n, sd = 1/sqrt(prec.y))

formula = y ~ f(w, model="meb",
        hyper = list(prec = list(param = c(1, 0.01))))

r = inla(formula, data = data.frame(y, w, s), 
        family = "gaussian")
