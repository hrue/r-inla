## When issue 14/15 is fixed, this file should run without error.

rate = log(1 / 10) # log of lambda
n = 200
n2 = 2*n -1 # if one removes -1 the code works
x <-  cumsum(rnorm(n2) * 0.05)
i_x = 1:n2
Q = toeplitz(c(2,-1, rep(0,(n2)-3),-1))
Q[1,1] = 1.01 #adding small 0 prior on first observation
Q[n2,n2] = 1
Q[n2,1] = Q[1,n2] = 0

A = matrix(data = 0, ncol = n2, nrow = n)
count = 1
for(i in 1:n){
	A[i,i] = 1
	A[i,count] = 1
	count = count + 2
}
rate_x = exp(A%*%x + rate)
time_1 <- rexp(n, rate = rate_x) # time of death
time_2 <- rexp(n, rate = exp(rate)) # censor time
event = (time_1 <= time_2 ) * 1
time = apply(cbind(time_1,time_2),1,min)
data = list(time = time, event = event, x = i_x)
formula <- inla.surv(time,event)  ~ f(x, model="generic", Cmatrix = Q)  - 1
model   <- inla(formula, family = "exponential" ,data = data, control.predictor= list(A=A), debug =FALSE )

formula <- inla.surv(time,event) ~ f(x, model="generic", Cmatrix = Q)  - 1
stack=
    inla.stack(data=list(time=time,event=event),
               A=A,
               effects=list(x=i_x))
inla.data = inla.stack.data(stack, Q=Q)
inla.A = inla.stack.A(stack)
model   <- inla(formula, family = "exponential" ,data = inla.data, control.predictor= list(A=inla.A), debug =FALSE, verbose=TRUE )

