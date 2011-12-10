## When issue 14/15 is fixed, this file should run without error.

require(INLA)
rate = log(1 / 10) # log of lambda
n = 200
x <-  cumsum(rnorm(n/2) * 0.05)
i_x = 1:(n/2)
Q = toeplitz(c(2,-1, rep(0,(n/2)-3),-1))
Q[1,1] = 1.1 #adding small 0 prior on first observation
Q[n/2,n/2] = 1
Q[n/2,1] = Q[1,n/2] = 0

## (both matrix and Matrix should work, but Matrix will be sparse
## automatically)
A = Matrix(0, ncol = n/2, nrow = n)
count = 1
for(i in 1:(n/2)){
    A[count,i] = 1
    count = count + 1
    A[count,i] = 1
    count = count + 1
}
rate_x = exp(as.vector(A%*%x) + rate)
time_1 <- rexp(n, rate = rate_x) # time of death
time_2 <- rexp(n, rate = exp(rate)) # censor time
event = (time_1 <= time_2 ) * 1
time = apply(cbind(time_1,time_2),1,min)
data = list(time = time, event = event, x = i_x)
formula <- inla.surv(time,event)  ~ f(x, model="generic", Cmatrix = Q)  - 1
model   <- inla(formula, family = "exponential" ,data = data, control.predictor= list(A=A), debug =FALSE, verbose=TRUE )

stack=
    inla.stack(data=list(time=time,event=event),
               A=A,
               effects=list(x=i_x),compress=FALSE)
inla.data = inla.stack.data(stack, Q=Q)
inla.A = inla.stack.data(stack)
model   <- inla(formula, family = "exponential" ,data = inla.data, control.predictor= list(A=inla.A), debug =FALSE, verbose=TRUE )
