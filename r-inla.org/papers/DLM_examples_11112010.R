
###########################
####   A toy example   ####
###########################

# Simulating observations and states from a first order DLM
# (based on code by Hedibert Lopes <http://faculty.chicagobooth.edu/hedibert.lopes/>
# ---------------------------------------------------------
set.seed(123456)

W = 0.5
V = 1.0
n = 100
x0 = 0.00
w  = rnorm(n,0,sqrt(W))
v  = rnorm(n,0,sqrt(V))
x  = y = rep(0,n)
x[1] = x0   + w[1]
y[1] = x[1] + v[1]
for (t in 2:n){
  x[t] = x[t-1] + w[t]
  y[t] = x[t]   + v[t]
}

m = n-1

# Plotting the simulated series
# -----------------------------
plot(y,type="l",xlab="time",ylab="")
lines(x,col=2)
legend("topright",legend=c("y","x"), col=c("black", "red"),lty=c(1,1),bty="n")
title(paste("V=",V," ; W=",W,sep=""))


# building the augmented model
# ----------------------------
Y <- matrix(NA, n+m, 2)
Y[1:n,     1] <- y              # actual observations
Y[1:m + n, 2] <- 0              # faked observations


# indices for the INLA library
# ----------------------------
i  <- c(1:n, 2:n)              # indices for x_t
j  <- c(rep(NA,n), 2:n -1)     # indices for x_{t-1}
w1 <- c(rep(NA,n), rep(-1,m))  # weights for j
l  <- c(rep(NA,n), 1:m)        # indices for w_t



# formulating the model
# ---------------------
formula <- Y ~ f(i, model="iid", initial=-10, fixed=T) +
               f(j, w1, copy="i") +
               f(l, model ="iid") -1


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r <- inla(formula, data = data.frame(i,j,w1,l),
          family = rep("gaussian", 2),
          control.family = list(list(), list(10, fixed=T)),
          control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

# elapsed time
r$cpu.used


## plotting the results
# ---------------------

# graph for observations (y)
rang <- range(r[[17]][1:n, 3:5], y)
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="y",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(y[1:n])
legend("topright", legend=c("simulated y_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")

# graph for states (x)
rang <- range(r[[12]][[1]][1:n, 4:6], x)
plot(r[[12]][[1]][1:n,2], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="X_t",xlab="time")
lines(r[[12]][[1]][1:n,4], col="blue", lty=3)
lines(r[[12]][[1]][1:n,6], col="blue", lty=3)
lines(x[1:n])
legend("topright", legend=c("simulated X_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")


# summary of posterior hyperparameters
# -------------------------------------
r$summary.hyperpar


# posterior density for precision parameters
# ------------------------------------------
par(mfrow=c(2,1))
plot(r[[23]][[1]], type="l", xlab="1/V", main="precision for Gaussian observations", ylab="")  # precision of V
  abline(v=1/V,col=2, lwd=5)
plot(r[[23]][[2]], type="l", xlab="1/W", main="precision for w_t", ylab="")  # precision of W
  abline(v=1/W,col=2, lwd=5)


# an alternative way to fit the same model using and RW1 process
i = 1:n
formula1 = y ~ f(i, model="rw1", constr=F) -1
r = inla(formula1, data = data.frame(i,y),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))


# graph for observations (y)
rang <- range(r[[17]][1:n, 3:5], y)
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="y",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(y[1:n])
legend("topright", legend=c("simulated y_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")

# graph for states (x)
rang <- range(r[[12]][[1]][1:n, 4:6], x)
plot(r[[12]][[1]][1:n,2], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="X_t",xlab="time")
lines(r[[12]][[1]][1:n,4], col="blue", lty=3)
lines(r[[12]][[1]][1:n,6], col="blue", lty=3)
lines(x[1:n])
legend("topright", legend=c("simulated X_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")


############################
####     EXAMPLE 1      ####
############################

# Simulating a seasonal DLM (seasonal length=12)
# ----------------------------------------------
set.seed(123456)

W1 = 0.001
V = 0.001
n = 100
a0 = 0.001
w1  = rnorm(n,0,sqrt(W1))
v  = rnorm(n,0,sqrt(V))
y = a  = rep(0,n)
a[1] = a0   + w1[1]
y[1] = a[1] + v[1]
for (t in 2:n){
  a[t] = a[t-1] + w1[t]
  y[t] = a[t]*cos(pi*(t-1)/6) + v[t]
}
y3 = y


# Plotting the seasonal DLM
# -------------------------
plot(y3,type="l",xlab="time",ylab="")
#lines(a1,col=2)
legend(80,9,legend=c("y","a_t"), col=c("black", "red"),lty=c(1,1),bty="n")
title(paste("V=",V," ; W1=",W1,sep=""))
#lines(a, col=2)

t = 1:n
cosw = cos(pi*(t-1)/6)


# fitting the model using and RW1 process
# ---------------------------------------
i=1:n
formula1 = y3 ~ f(i, cosw, model="rw1", constr=F) -1
r = inla(formula1, data = data.frame(i,y3),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

# graph for observations (y)
rang <- range(r[[17]][1:n, 3:5], y3)
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="y",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(y3[1:n])
legend("topright", legend=c("simulated y_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")

# graph for states (x)
rang <- range(r[[12]][[1]][1:n, 4:6], a)
plot(r[[12]][[1]][1:n,2], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="X_t",xlab="time")
lines(r[[12]][[1]][1:n,4], col="blue", lty=3)
lines(r[[12]][[1]][1:n,6], col="blue", lty=3)
lines(a[1:n])
legend("topright", legend=c("simulated X_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")


# hyperparameter estimation values
# ---------------------------------
r$summary.hyperpar


# posterior density for precision parameters
# ------------------------------------------
par(mfrow=c(2,2))
plot(r[[23]][[1]], type="l", xlab="1/V", main="precision for y", ylab="")  # precision of V
  abline(v=1/V,col=2, lwd=3)
plot(r[[23]][[2]], type="l", xlab="1/W1", main="precision for x_1", ylab="")  # precision of W
  abline(v=1/W1,col=2, lwd=3)



############################
####     EXAMPLE 2      ####
############################

## Simulating a Poisson dynamic regression
## (dynamic intercept and two dynamic regrressors)
## -----------------------------------------------
set.seed(123456)

n <- 300
x1 <- runif(n)
x2 <- runif(n)
W0 = 0.005
W1 = 0.01
W2 <- 0.005
w0 <- rnorm(n,0,sqrt(W0))
w1 <- rnorm(n,0,sqrt(W1))
w2 <- rnorm(n,0,sqrt(W2))
y <- b0 <- b1 <- b2 <- numeric(n)
b0[1] <- 3 + w0[1]
b1[1] <- 0 + w1[1]
b2[1] <- 0 + w2[1]
for (t in 2:n){ 
  b0[t] <- b0[t-1] + w0[t]
  b1[t] <- b1[t-1] + w1[t]
  b2[t] <- b2[t-1] + w2[t]
}
y <- rpois(n,exp(b0 + b1*x1 + b2*x2))

## Plotting the dynamic regression model
## -------------------------------------
par(mfrow=c(2,1))
plot.ts(y,xlab="time",ylab="",ylim=range(c(y,b1)))
lines(exp(b0 + b1*x1 + b2*x2),col=2)
legend("topleft",legend=c("y","lambda"),
       col=c("black", "red"),lty=1,bty="n")
plot.ts(b0,xlab="time",ylab="",ylim=range(c(b0,b1,b2)))
lines(b1,col=4)
lines(b2,col=2)
legend("bottomleft",legend=c("b0","b1","b2"),
       col=c("black", "blue","red"),lty=1,bty="n")
title(paste("W1=",W1, ", W2=", W2,sep=""))


# defining indices for rw1 coefficients
# -------------------------------------
id <- id1 <- id2 <- 1:n


# formulating the model
# ---------------------
formula <- y ~ f(id, model="rw1",initial=5,constr=F) + 
               f(id1, x1, model="rw1",initial=5,constr=F) + 
               f(id2, x2, model="rw1",initial=5,constr=F) -1

formula2 <- y ~ f(id, model="rw1",param=c(2,0.01),initial=5,constr=F) + 
               f(id1, x1, model="rw1",param=c(2,0.02),initial=5,constr=F) + 
               f(id2, x2, model="rw1",param=c(2,0.01),initial=5,constr=F) -1

formula3 <- y ~ f(id, model="rw1",param=c(0.1,0.0005),initial=5,constr=F) + 
               f(id1, x1, model="rw1",param=c(0.1,0.001),initial=5,constr=F) + 
               f(id2, x2, model="rw1",param=c(0.1,0.0005),initial=5,constr=F) -1

require(INLA)
r = inla(formula, family="poisson", data = data.frame(id,id1,id2,y),
         control.predictor=list(compute=TRUE))

r2 = inla(formula2, family="poisson", data = data.frame(id,id1,id2,y),
         control.predictor=list(compute=TRUE))

r3 = inla(formula3, family="poisson", data = data.frame(id,id1,id2,y),
         control.predictor=list(compute=TRUE))


## plotting the results

## graph for y
par(mfrow=c(2,2))
rang <- range(r[[17]][1:n, 3:5], y)
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="y",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(y)
legend("topleft", legend=c("simulated","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("dynamic regression model")

## graph for states (b0_t)
rang <- range(r[[12]][[1]][1:n,4:6], b0[1:n])
plot(r[[12]][[1]][1:n,2], type="l", 
     ylim=rang, col="red", xlim=c(1,n),ylab="b0",xlab="time")
lines(r[[12]][[1]][1:n,4], col="blue", lty=3)
lines(r[[12]][[1]][1:n,6], col="blue", lty=3)
lines(b0[1:n])
legend("topleft", legend=c("simulated","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
##title("a toy example")

## graph for states (b1_t)
rang <- range(r[[12]][[2]][1:n,4:6], b1[1:n])
plot(r[[12]][[2]][1:n,2], type="l", 
     ylim=rang, col="red", xlim=c(1,n),ylab="b1",xlab="time")
lines(r[[12]][[2]][1:n,4], col="blue", lty=3)
lines(r[[12]][[2]][1:n,6], col="blue", lty=3)
lines(b1[1:n])
legend("bottomleft", legend=c("simulated","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
##title("a toy example")

## graph for states (b2_t)
rang <- range(r[[12]][[3]][1:n,4:6], b2[1:n])
plot(r[[12]][[3]][1:n,2], type="l", 
     ylim=rang, col="red", xlim=c(1,n),ylab="b2",xlab="time")
lines(r[[12]][[3]][1:n,4], col="blue", lty=3)
lines(r[[12]][[3]][1:n,6], col="blue", lty=3)
lines(b2[1:n])
legend("topleft", legend=c("simulated","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
##title("a toy example")


# hyperparamerter estimation values
# ---------------------------------
r$summary.hyperpar
r2$summary.hyperpar
r3$summary.hyperpar


# posterior densities for precision parameters
# --------------------------------------------
par(mfrow=c(2,2))
rang1 <- range(r[[23]][[1]][,2],r2[[23]][[1]][,2],r3[[23]][[1]][,2])
plot(r[[23]][[1]], type="l", ylab="",xlab="1/W0" , main="precision for w_0", ylim=rang1)  # precision of W0
  lines(r2[[23]][[1]],lty=3)
  lines(r3[[23]][[1]],lty=2)
abline(v=1/W0,col=2, lwd=3)
rang2 <- range(r[[23]][[2]][,2],r2[[23]][[2]][,2],r3[[23]][[2]][,2])
plot(r[[23]][[2]], type="l", main="precision for w_1", ylab="" ,ylim=rang2, xlab="1/W1")  
  lines(r2[[23]][[2]],lty=3)
  lines(r3[[23]][[2]],lty=2)
abline(v=1/W1,col=2, lwd=3)
rang3 <- range(r[[23]][[3]][,2],r2[[23]][[3]][,2],r3[[23]][[3]][,2])
plot(r[[23]][[3]], type="l", main="precision for w_2", ylab="" ,ylim=rang3, xlab="1/W2")  # precision of W2
  lines(r2[[23]][[3]],lty=3)
  lines(r3[[23]][[3]],lty=2)
legend("topright", legend=c("with informative prior","with vague prior","with default prior"), col=c("black", "black", "black"),lty=c(3,2,1),bty="n",cex=1)
abline(v=1/W2,col=2, lwd=3)



############################
####     EXAMPLE 3      ####
############################

# Simulating observations and states from a second order DLM 
# (based on code by Hedibert Lopes <http://faculty.chicagobooth.edu/hedibert.lopes/>
# ----------------------------------------------------------
set.seed(1)

W1 = 0.0001
W2 = 0.0001
V1 = 0.01
n  = 100
x0 = c(0.00,0.00)
w  = cbind(rnorm(n,0,sqrt(W1)),rnorm(n,0,sqrt(W2)))
v  = rnorm(n,0,sqrt(V1))
y = rep(0,n)
x = matrix(0,n,2)
x[1,2] = x0[2]          + w[1,2]
x[1,1] = x0[1] - x0[2]  + w[1,1]
y[1]   = x[1,1]         + v[1]
for (t in 2:n){
  x[t,2] = x[t-1,2]            + w[t,2]
  x[t,1] = x[t-1,1] - x[t-1,2] + w[t,1]
  y[t]   = x[t,1]              + v[t]
}

m=n-1


# Plotting simulated series
# -------------------------
plot(y,type="l",xlab="time",ylab="")
lines(x[,1],col=2)
lines(x[,2],col=4)
legend(0,13,legend=c("y","x1","x2"), col=c("black", "red","blue"),lty=c(1,1,1),bty="n")
title(paste("V=",V1," ; W1=",W1," ; W2=",W2,sep=""))


# building the augmented model
# ----------------------------
Y <- matrix(NA, n+2*m, 3)
Y[1:n,     1] <- y                  # actual observations (y)
Y[1:m + n, 2] <- 0                  # faked observations (x1)
Y[1:m + (n+m), 3] <- 0              # faked observations (x2)


# indices for the INLA library
# ----------------------------
i       = c(1:n, 2:n, rep(NA,m))             # x_{1,t}
j       = c(rep(NA,n), 2:n -1, rep(NA,m))    # x_{1,t-1}
weight1 = c(rep(NA,n), rep(-1,m), rep(NA,m)) # weights for j
l       = c(rep(NA,n+m), 2:n)                # x_{2,t}
k       = c(rep(NA,n), 2:n -1, 2:n -1)       # x_{2,t-1} for the first state equation
weight2 = c(rep(NA,n), rep(1,m), rep(-1,m))  # weights for k
q       = c(rep(NA,n), 1:m, rep(NA,m))       # w_{1,t}
s       = c(rep(NA,n+m), 1:m)                # w_{2,t}


# formulating the model
# ---------------------
formula = Y ~ f(i, model="iid", initial=-10, fixed=T) +
              f(j, weight1, copy="i") + 
              f(l, model="iid", initial=-10, fixed=T) + 
              f(k, weight2, copy="l") +
              f(q, model ="iid") + 
              f(s, model ="iid") -1


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r = inla(formula, data = data.frame(i,j,weight1,k,weight2,l,q,s),
         family = rep("gaussian", 3),
#          control.inla=list(strategy="laplace", int.strategy="grid"),
         control.family = list(list(),list(initial=10, fixed=T),list(initial=10, fixed=T)),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))


# elapsed time
r$cpu.used


## plotting the results
# ---------------------

# graph for observations (y)
rang <- range(r[[17]][1:n, 3:5], y)
plot(r[[17]][1:n,1], type="l", 
ylim=rang, col="red", xlim=c(1,n),ylab="y_t",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(y[1:n])
legend("topleft",legend=c("simulated y_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("second order polynomial model")

# graph for states (x1)
rang <- range(r[[12]][[1]][1:n,4:6], x[1:n,1])
plot(r[[12]][[1]][1:n,2], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="x_1t",xlab="time")
lines(r[[12]][[1]][1:n,4], col="blue", lty=3)
lines(r[[12]][[1]][1:n,6], col="blue", lty=3)
lines(x[1:n,1])
legend("topleft",legend=c("simulated X_1","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")

# graph for states (x2)
rang <- range(r[[12]][[2]][1:(n-1),4:6], x[2:n,2])
plot(r[[12]][[2]][1:(n-1),2], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="x_2t",xlab="time")
lines(r[[12]][[2]][1:(n-1),4], col="blue", lty=3)
lines(r[[12]][[2]][1:(n-1),6], col="blue", lty=3)
lines(x[2:n,2])
legend("bottomleft",legend=c("simulated X_2","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")


# hyperparameter estimation values
# ---------------------------------
r$summary.hyperpar


# posterior density for precision parameters
# ------------------------------------------
par(mfrow=c(2,2))
plot(r[[23]][[1]], type="l", xlab="1/V", main="precision for observations", ylab="")  # precision of V
  abline(v=1/V1,col=2, lwd=3)
plot(r[[23]][[2]], type="l", xlab="1/W1", main="precision for w_1t", ylab="")  # precision of W
  abline(v=1/W1,col=2, lwd=3)
plot(r[[23]][[3]], type="l", xlab="1/W2", main="precision for w_2t", ylab="")  # precision of W
  abline(v=1/W2,col=2, lwd=3)



############################
####     EXAMPLE 4      ####
############################

# Simulating a seasonal dynamic model with harmonics

set.seed(123456)

W1 = 0.002
W2 = 0.001
V = 0.01
n = 100
a0 = 0.1
b0 = 0.1
w1  = rnorm(n,0,sqrt(W1))
w2  = rnorm(n,0,sqrt(W2))
v  = rnorm(n,0,sqrt(V))
y = a  = b = rep(0,n)
a[1] = cos(pi/6)*a0 + sin(pi/6)*b0 + w1[1]
b[1] = -sin(pi/6)*a0 + cos(pi/6)*b0 + w2[1]
y[1] = a[1] + v[1]
for (t in 2:n){
  a[t] = cos(pi/6)*a[t-1] + sin(pi/6)*b[t-1] + w1[t]
  b[t] = -sin(pi/6)*a[t-1] + cos(pi/6)*b[t-1] + w2[t]
  y[t] = a[t] + v[t]
}


# Plotting the seasonal DLM
# -------------------------
plot(y,type="l",xlab="time",ylab="")
lines(a,col=2)
lines(b,col=4)
legend("topleft",legend=c("y","a_t","b_t"), col=c("black", "red","blue"),lty=c(1,1,1),bty="n")
title(paste("V=",V," ; W1=",W1," ; W2=",W2,sep=""))


m = n-1
cosw = cos(pi/6)
sinw = sin(pi/6)

# building the augmented model
# ----------------------------
Y <- matrix(NA, n+2*m, 3)
Y[1:n,         1] <- y              # actual observations (y)
Y[1:m + n,     2] <- 0              # faked observations (at)
Y[1:m + (n+m), 3] <- 0              # faked observations (bt)


# indices for the INLA library
# ----------------------------
i       = c(1:n, 2:n, rep(NA,m))                  # a_t
weight0 = c(rep(1,n+m), rep(NA,m))                # weights for i
j       = c(rep(NA,n), 2:n -1, 2:n -1)            # a_{t-1}
weight1 = c(rep(NA,n), rep(-cosw,m), rep(sinw,m)) # weights for j
l       = c(rep(NA,n+m), 2:n)                     # b_t
weight2 = c(rep(NA,n+m), rep(1,m))                # weights for l
o       = c(rep(NA,n), 2:n -1, 2:n -1)            # b_{t-1}
weight3 = c(rep(NA,n), rep(-sinw,m), rep(-cosw,m))# weights for o
q       = c(rep(NA,n), 2:n, rep(NA,m))            # w_{1,t}
s       = c(rep(NA,n+m), 2:n)                     # w_{2,t}


# formulating the model
# ---------------------
formula = Y ~ f(i,weight0, model="iid", initial=-10, fixed=T) +
              f(j, weight1, copy="i") +
              f(l, weight2, model="iid", initial=-10, fixed=T) +
              f(o, weight3, copy="l") +
              f(q, model ="iid") + f(s, model ="iid") -1


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r = inla(formula, data = data.frame(i,weight0,j,weight1,l,weight2,o,weight3,q,s),
         family = rep("gaussian", 3),
#          control.inla=list(strategy="laplace", int.strategy="grid"),
         control.family = list(list(param=c(1,0.01),initial=log(20)),list(initial=log(10), fixed=T),list(initial=log(20), fixed=T)),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

# elapsed time
r$cpu.used


## plotting the results
# ---------------------

# graph for observations (y)
rang <- range(r[[17]][1:n, 3:5], y)
plot(r[[17]][1:n,1], type="l", 
ylim=rang, col="red", xlim=c(1,n),ylab="y",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(y[1:n])
legend("topleft",legend=c("simulated y_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("second order polynomial model")

# graph for states (a_t)
rang <- range(r[[12]][[1]][1:n,4:6], a[1:n])
plot(r[[12]][[1]][1:n,2], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="a_t",xlab="time")
lines(r[[12]][[1]][1:n,4], col="blue", lty=3)
lines(r[[12]][[1]][1:n,6], col="blue", lty=3)
lines(a[1:n])
legend("topleft", legend=c("simulated a_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")

# graph for states (b_t)
rang <- range(r[[12]][[2]][1:m,c(4,6)], b[1:m])
plot(r[[12]][[2]][1:m,2], type="l", 
 ylim=rang, col="red", xlim=c(2,n),ylab="b_t",xlab="time")
lines(r[[12]][[2]][1:m,4], col="blue", lty=3)
lines(r[[12]][[2]][1:m,6], col="blue", lty=3)
lines(b[1:m])
legend("bottomleft", legend=c("simulated b_t","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("a toy example")


# hyperparameter estimation values
# ---------------------------------
r$summary.hyperpar


# posterior densities for precision parameters
# --------------------------------------------
par(mfrow=c(2,2))
plot(r[[23]][[1]], type="l", xlab="1/V", main="precision for y", ylab="")  # precision of V
  abline(v=1/V,col=2, lwd=3)
plot(r[[23]][[2]], type="l", xlab="1/W1", main="precision for a_t", ylab="")  # precision of W
  abline(v=1/W1,col=2, lwd=3)
plot(r[[23]][[3]], type="l", xlab="1/W2", main="precision for b_t", ylab="")  # precision of W
  abline(v=1/W2,col=2, lwd=3)



############################
####     EXAMPLE 5      ####
############################

# Here we simulated data from a non-stationary first-order Gaussian dynamic 
# spatio-temporal model with one spatially structured covariate 
# (Vivar & Ferreira, 2009) defined as follows:
#
# Observational equation:
#  y_{i,t} = x_{i,t} + \beta1*Z1_{i,t} + w1_{i,t},    i=1,...,n,  t=1,...,k
# System equations:
#  x_{i,t} = x_{i,t-1} + w2_{i,t},                    i=1,...,n,  t=2,...,k
#  Z1_{i,t} = Z1_{i,t-1} + w3_{i,t},                  i=1,...,n,  t=2,...,k
#
# We can re-write the system equation model as
#  0 = x_{i,t} - x_{i,t-1} - w2_{i,t},                i=1,...,n,  t=2,...,k
#  0 = Z1_{i,t} - Z1_{i,t-1} - w3_{i,t},              i=1,...,n,  t=2,...,k
# where 
#  w1_{,t} \sim N(0, \Sigma_1)
#  w2_{,t} \sim N(0, \Sigma_2)
#  w3_{,t} \sim N(0, \Sigma_3)
#
# we use \Sigma_j^{-1} = \tau_j(I - \phi_j C), j = {1, 2}
# where
#  D_j =  diag\{d_1, d_2, ..., d_n\}
#  d_i    is the number of neighbours of area i
#  \tau_j is the precision parameter
#  \phi_j is the spatial correlation parameter
#  C      is the standardized adjacency matrix (here each row sums one)


### simulating a data set from this model

## Loading map of Eire (it has 26 areas)
require(spdep)
ncfile <- system.file("etc/shapes/eire.shp", package="spdep")[1]
nc <- readShapePoly(ncfile)

# building the structure matrix (C)
nc.nb <- poly2nb(nc)                      
d <- sapply(nc.nb, length)                   # vector with number of neighbors
C <- diag(d) - nb2mat(nc.nb, style="B")      # structure matrix

n <- length(d)

# simulated values for tau_i and phi_i (i=1,2)
tau <- c(30, 50, 50)
phi <- c(0.8, 0.9, 0.8)  

# building the precision matrix
lamb.max <- max(eigen(C, only.values=TRUE)$values) # maximum eigenvalue of C matrix
Q1 <- (diag(n)-phi[1]/lamb.max*C)
Q2 <- (diag(n)-phi[2]/lamb.max*C)
Q3 <- (diag(n)-phi[3]/lamb.max*C)

myrmvnorm <- function(n, mu, S) 
  sweep(matrix(rnorm(n*nrow(S)), n)%*%chol(S), 2, mu)

# defining the length of time series (number of years)
k <- 100

set.seed(1)

# simulating obsevational and innovation errors
w1 <- t(myrmvnorm(k, rep(0,n), solve(tau[1]*Q1)))
w2 <- t(myrmvnorm(k, rep(0,n), solve(tau[2]*Q2)))
w3 <- t(myrmvnorm(k, rep(0,n), solve(tau[3]*Q2)))

# generating the time series for observations and states
beta1 <- matrix(rep(runif(n),k), n, k, byrow=TRUE)
y <- x <- Z <- matrix(0, n, k)
x[,1] <- w2[,1]
Z[,1] <- w3[,1]
for (i in 2:k){
  x[,i] <- x[,i-1] + w2[,i]
  Z[,i] <- Z[,i-1] + w3[,i]
}
y <- x + beta1*Z + w1

### defining the Cmatrix to use with model='generic1' for v and w
st.cmat <- kronecker(C, diag(k))
c.mat <- list(i=unlist(apply(st.cmat!=0, 1, which)),
              j=rep(1:nrow(st.cmat), rowSums(st.cmat!=0)), 
              values=st.cmat[st.cmat!=0])
##rm(st.cmat)


### building the augmented model
### ----------------------------
nd <- n*k
Y <- matrix(NA, nd*3-2*n, 3)
Y[1:nd         , 1] <- as.vector(t(y))
Y[1:(nd-n) + nd, 2] <- 0
Y[1:(nd-n) + 2*nd-n, 3] <- 0


### indices for the f() function
### ----------------------------
id1   <- (1:nd)[-((1:n)*k)]
id2   <- (1:nd)[-c(1,((1:(n-1))*k)+1)]
ix    <- c(1:nd, id2, rep(NA,nd-n))                        ## x_{,t}
ixb   <- c(rep(NA,nd), id1, rep(NA,nd-n))                  ## x_{,t-1}
w.ixb <- c(rep(NA,nd), rep(-1,nd-n), rep(NA,nd-n))         ## weights for x_{,t-1}
iZ    <- c(1:nd, rep(NA,nd-n), id2)                        ## Z_{,t}
iZb   <- c(rep(NA,nd), rep(NA,nd-n), id1)                  ## Z_{,t-1}
w.iZb <- c(rep(NA,nd), rep(NA,nd-n), rep(-1,nd-n))         ## weights for Z_{,t-1}
iw1   <- c(1:nd, rep(NA,2*(nd-n)))                         ## w1_{,t} 
iw2   <- c(rep(NA,nd), id2, rep(NA,nd-n))                  ## w2_{,t}
iw3   <- c(rep(NA,nd), rep(NA,nd-n), id2)                  ## w3_{,t}
beta  <- c(as.vector(t(beta1)), rep(NA,nd-n), rep(1,nd-n)) ## covariate
#w.iw2 <- c(rep(NA,nd), rep(-1,nd-n), rep(NA,nd-n))
#w.iw3 <- c(rep(NA,nd), rep(NA,nd-n), rep(-1,nd-n))


## formulating the model
## ---------------------
# with default prior for \tau and \phi
formula = Y ~ f(iZ, beta, model="iid", initial=-10, fixed=T) +
              f(iZb, w.iZb, copy="iZ") +
              f(iw1, model="generic1", Cmatrix=c.mat) +
              f(ix, model="iid", initial=-10, fixed=T) +
              f(ixb, w.ixb, copy="ix") + 
              f(iw2, model="generic1", Cmatrix=c.mat) + 
              f(iw3, model="generic1", Cmatrix=c.mat) -1 

# with informative prior for \tau and default prior for \phi
formula1 = Y ~ f(iZ, beta, model="iid", initial=-10, fixed=T) +
               f(iZb, w.iZb, copy="iZ") +
               f(iw1, model="generic1", Cmatrix=c.mat, param=c(1,1/30, 0,0.001)) +
               f(ix, model="iid", initial=-10, fixed=T) +
               f(ixb, w.ixb, copy="ix") + 
               f(iw2, model="generic1", Cmatrix=c.mat, param=c(1,1/50, 0,0.001)) + 
               f(iw3, model="generic1", Cmatrix=c.mat, param=c(1,1/50, 0,0.001)) -1 

# with vague prior for \tau and default prior for \phi
formula2 = Y ~ f(iZ, beta, model="iid", initial=-10, fixed=T) +
               f(iZb, w.iZb, copy="iZ") +
               f(iw1, model="generic1", Cmatrix=c.mat, param=c(4,0.13333, 0,0.001)) +
               f(ix, model="iid", initial=-10, fixed=T) +
               f(ixb, w.ixb, copy="ix") + 
               f(iw2, model="generic1", Cmatrix=c.mat, param=c(4,0.08, 0,0.001)) + 
               f(iw3, model="generic1", Cmatrix=c.mat, param=c(4,0.08, 0,0.001)) -1 



# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r = inla(formula, data = data.frame(beta,ix,ixb,w.ixb,iw1,iw2,iw3,iZ,iZb,w.iZb),
         family = rep("gaussian",3), 
         control.family = list(list(initial=10, fixed=T),list(initial=10, fixed=T),list(initial=10, fixed=T)),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

r1 = inla(formula1, data = data.frame(beta,ix,ixb,w.ixb,iw1,iw2,iw3,iZ,iZb,w.iZb),
         family = rep("gaussian",3), 
         control.family = list(list(initial=10, fixed=T),list(initial=10, fixed=T),list(initial=10, fixed=T)),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

r2 = inla(formula2, data = data.frame(beta,ix,ixb,w.ixb,iw1,iw2,iw3,iZ,iZb,w.iZb),
         family = rep("gaussian",3), 
         control.family = list(list(initial=10, fixed=T),list(initial=10, fixed=T),list(initial=10, fixed=T)),
         control.predictor=list(compute=TRUE, cdf=c(.025, .975)))


# elapsed time
r$cpu.used
r1$cpu.used
r2$cpu.used

# graph for observations (y)
rang <- range(r[[17]][1:nd, 3:5], y)
plot.ts(r[[17]][1:nd,1], type="l", ylim=rang, col="red", ylab="y",xlab="time")
lines(r[[17]][1:nd,3], col="blue", lty=3)
lines(r[[17]][1:nd,5], col="blue", lty=3)
lines(as.vector(t(y)))
legend("topleft",legend=c("observed","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("second order polynomial model")

#names(r[[12]])

par(mfrow=c(2,2))
plot(as.vector(t(w1)), r[[12]]$iw1$mean)
plot(as.vector(t(w2)), r[[12]]$iw2$mean)
plot(as.vector(t(w3)), r[[12]]$iw3$mean)
par(mfrow=c(2,2))
plot(as.vector(t(x)), r[[12]]$ix$mean)
plot(as.vector(t(Z)), r[[12]]$iZ$mean)


## graph for observations in the first 4 areas (i=1,2,3,4)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
for (i in 1:4) {
  id <- 1:k + (i-1)*k
  rang <- range(r[[17]][id,3:5], y[i,])
  plot.ts(r[[17]][id,1], ylim=rang, col="red", ylab="y_it",xlab="time")
  lines(r[[17]][id,3], col="blue", lty=3)
  lines(r[[17]][id,5], col="blue", lty=3)
  lines(y[i,])
  legend("topleft", legend=c("simulated y","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
}

## graph for states (X_{i,t}) in the first 4 areas (i=1,2,3,4)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
for (i in 1:4) {
  id <- 1:k + (i-1)*k
  rang <- range(r[[12]]$ix[id, 4:6], x[i,])
  plot.ts(r[[12]]$ix$mean[id], ylim=rang, col="red", ylab="X_it",xlab="time")
  lines(r[[12]]$ix[id,4], col="blue", lty=3)
  lines(r[[12]]$ix[id,6], col="blue", lty=3)
  lines(x[i,])
  legend("topleft", legend=c("simulated X","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
}

## graph for states (Z_{i,t}) in the first 4 areas (i=1,2,3,4)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
for (i in 1:4) {
  id <- 1:k + (i-1)*k
  rang <- range(r[[12]]$iZ[id, 4:6], Z[i,])
  plot.ts(r[[12]]$iZ$mean[id], ylim=rang, col="red", ylab="Z_it",xlab="time")
  lines(r[[12]]$iZ[id,4], col="blue", lty=3)
  lines(r[[12]]$iZ[id,6], col="blue", lty=3)
  lines(Z[i,])
  legend("topleft", legend=c("simulated Z","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
}


## graph of observations y_t for a given area and its neighbors
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
area=18
for (i in c(area, nc.nb[[area]])) {
  id <- 1:k + (i-1)*k
  rang <- range(r[[17]][id,3:5], y[i,])
  plot.ts(r[[17]][id,1], ylim=rang, col="red", ylab="y_t",xlab="time")
  lines(r[[17]][id,3], col="blue", lty=3)
  lines(r[[17]][id,5], col="blue", lty=3)
  lines(y[i,])
}
legend("bottomleft",legend=c("simulated","posterior mean", "95% CI"),
       col=c("black","red","blue"),lty=c(1,1,2),bty="n")


## graph of states X_t for a given area and its neighbors
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
area=18
for (i in c(area, nc.nb[[area]])) {
  id <- 1:k + (i-1)*k
  rang <- range(r[[12]]$ix[id, 4:6], x[i,])
  plot.ts(r[[12]]$ix$mean[id], ylim=rang, col="red", ylab="X_t",xlab="time")
  lines(r[[12]]$ix[id,4], col="blue", lty=3)
  lines(r[[12]]$ix[id,6], col="blue", lty=3)
  lines(x[i,])
}
legend("bottomleft",legend=c("simulated","posterior mean", "95% CI"),
       col=c("black","red","blue"),lty=c(1,1,2),bty="n")

## graph of states Z_t for a given area and its neighbors
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
area=18
for (i in c(area, nc.nb[[area]])) {
  id <- 1:k + (i-1)*k
  rang <- range(r[[12]]$iZ[id, 4:6], Z[i,])
  plot.ts(r[[12]]$iZ$mean[id], ylim=rang, col="red", ylab="Z_t",xlab="time")
  lines(r[[12]]$iZ[id,4], col="blue", lty=3)
  lines(r[[12]]$iZ[id,6], col="blue", lty=3)
  lines(Z[i,])
}
legend("bottomleft",legend=c("simulated","posterior mean", "95% CI"),
       col=c("black","red","blue"),lty=c(1,1,2),bty="n")


# hyperparameter estimation values
# ---------------------------------
r$summary.hyperpar
r1$summary.hyperpar
r2$summary.hyperpar


# posterior densities for precision parameters
# --------------------------------------------
par(mfrow=c(2,2))
rang1 <- range(r[[23]][[1]][,2],r1[[23]][[1]][,2],r2[[23]][[1]][,2])
plot(r[[23]][[1]], type="l", main=expression(tau[W1]), ylab="",xlab="" ,ylim=rang1)  # precision of V
  lines(r1[[23]][[1]],lty=3)
  lines(r2[[23]][[1]],lty=2)
abline(v=tau[1],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang3 <- range(r[[23]][[3]][,2],r1[[23]][[3]][,2],r2[[23]][[3]][,2])
plot(r[[23]][[3]], type="l", main=expression(tau[W2]), ylab="" ,ylim=rang3,xlab="")  # precision of V
  lines(r1[[23]][[3]],lty=3)
  lines(r2[[23]][[3]],lty=2)
abline(v=tau[2],col=2, lwd=3)
legend("topright", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang5 <- range(r[[23]][[5]][,2],r1[[23]][[5]][,2],r2[[23]][[5]][,2])
plot(r[[23]][[5]], type="l", main=expression(tau[W3]), ylab="" ,ylim=rang5,xlab="")  # precision of V
  lines(r1[[23]][[5]],lty=3)
  lines(r2[[23]][[5]],lty=2)
abline(v=tau[3],col=2, lwd=3)
legend("topright", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)


par(mfrow=c(2,2))
rang2 <- range(r[[23]][[2]][,2],r1[[23]][[2]][,2],r2[[23]][[2]][,2])
plot(r[[23]][[2]], type="l", main=expression(phi[W1]), ylab="" ,ylim=rang2,xlab="")  
  lines(r1[[23]][[2]],lty=3)
  lines(r2[[23]][[2]],lty=2)
abline(v=phi[1],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang4 <- range(r[[23]][[4]][,2],r1[[23]][[4]][,2],r2[[23]][[4]][,2])
plot(r[[23]][[4]], type="l", main=expression(phi[W2]), ylab="" ,ylim=rang4,xlab="")  
  lines(r1[[23]][[4]],lty=3)
  lines(r2[[23]][[4]],lty=2)
abline(v=phi[2],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang6 <- range(r[[23]][[6]][,2],r1[[23]][[6]][,2],r2[[23]][[6]][,2])
plot(r[[23]][[6]], type="l", main=expression(phi[W3]), ylab="" ,ylim=rang6,xlab="")  
  lines(r1[[23]][[6]],lty=3)
  lines(r2[[23]][[6]],lty=2)
abline(v=phi[3],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)


require(RColorBrewer)
c4 <- c(rev(brewer.pal(3, "Blues"))[1:2], brewer.pal(3, "Reds")[2:3])

### map of y at times 34, 64 and 91
postscript("D:\\meus artigos\\INLA-DLMs\\figures/ex5Z_mapy_346491.eps", 
  horizontal=FALSE, width=17, height=5)
par(mfrow=c(3,2), mar=c(0,0,0,0))
ind1 = c(34,64,91)
for (t in ind1) {
  q <- quantile(c(y[,t], r[[17]][seq(1, nd, k)+t-1,1]), 0:4/4)
  plot(nc, col=c4[findInterval(y[,t], q, T)])
#title(t, "areas")
  plot(nc, col=c4[findInterval(r[[17]][seq(1, nd, k)+t-1,1], q, T)])
#title("areas")
  legend("topleft", leglabs(format(q,dig=1),bet="to"), fill=c4, bty="n", ncol=1)
}
dev.off()

### map of x1 at times 34, 64 and 91
postscript("D:\\meus artigos\\INLA-DLMs\\figures/ex5Z_mapx1_346491.eps", 
  horizontal=FALSE, width=17, height=5)
par(mfrow=c(3,2), mar=c(0,0,0,0))
ind1 = c(34,64,91)
for (t in ind1) {
  q <- quantile(c(x[,t], r[[12]]$ix$mean[seq(1,nd,k)+t-1]), 0:4/4)
  plot(nc, col=c4[findInterval(x[,t], q, T)])
  plot(nc, col=c4[findInterval(r[[12]]$ix$mean[seq(1, nd, k)+t-1], q, T)])
  legend("topleft", leglabs(format(q,dig=1),bet="to"), fill=c4, bty="n", ncol=1)
}
dev.off()


### map of Z at times 34, 64 and 91
postscript("D:\\meus artigos\\INLA-DLMs\\figures/ex5Z_mapZ_346491.eps", 
  horizontal=FALSE, width=17, height=5)
par(mfrow=c(3,2), mar=c(0,0,0,0))
ind1 = c(34,64,91)
for (t in ind1) {
  q <- quantile(c(Z[,t], r[[12]]$iZ$mean[seq(1,nd,k)+t-1]), 0:4/4)
  plot(nc, col=c4[findInterval(Z[,t], q, T)])
  plot(nc, col=c4[findInterval(r[[12]]$iZ$mean[seq(1, nd, k)+t-1], q, T)])
  legend("topleft", leglabs(format(q,dig=1),bet="to"), fill=c4, bty="n", ncol=1)
}
dev.off()



############################
####     EXAMPLE 6      ####
############################

# Here we simulate data from a non-stationary second-order Gaussian dynamic 
# spatio-temporal model (Vivar & Ferreira, 2009) defined as follows:
#
# Observational equation:
#  y_{i,t} = x1_{i,t} + w1_{i,t},                     i=1,...,n, t=1,...,T
# System equations:
#  x1_{i,t} = x1_{i,t-1} + x2_{i,t-1} + w2_{i,t},     i=1,...,n, t=2,...,T
#  x2_{i,t} = x2_{i,t-1} + w3_{i,t},                  i=1,...,n, t=2,...,T
#
# We can re-write the system equations as
#  0 = x1_{i,t} - x1_{i,t-1} - x2_{i,t-1} - w2_{i,t}, i=1,...,n, t=2,...,T
#  0 = x2_{i,t} - x2_{i,t-1} - w3_{i,t},              i=1,...,n, t=2,...,T
#
# where 
#  w1_{,t} \sim N_n(0, \Sigma_1)
#  w2_{,t} \sim N_n(0, \Sigma_2)
#  w3_{,t} \sim N_n(0, \Sigma_3)
#
# we use \Sigma_j^{-1} = D_j(I - \phi_j C), j = {1, 2, 3}
# where
#  D_j =  \tau_j diag\{d_1, d_2, ..., d_n\}
#  d_i    is the number of neighbours of area i
#  \tau_j is the precision parameter
#  \phi_j is the spatial correlation parameter
#  C      is the standardized adjacency matrix (here each row sums one)

## simulating the data set

# Loading map of North Carolina (it has 100 areas)
require(spdep)
ncfile <- system.file("etc/shapes/sids.shp", package="spdep")[1]
nc <- readShapePoly(ncfile)

# building the structure matrix (C)
nc.nb <- poly2nb(nc)                      
d <- sapply(nc.nb, length)                   # vector with number of neighbors
C <- diag(d) - nb2mat(nc.nb, style="B")      # structure matrix

n <- length(d)

# simulated values for tau_i and phi_i (i=1,2,3)
tau <- c(30, 50, 70)
phi <- c(0.8, 0.9, 0.9)

# building the precision matrix
lamb.max <- max(eigen(C, only.values=TRUE)$values) # maximum eigenvalue of C matrix
Q1 <- (diag(n)-phi[1]/lamb.max*C)
Q2 <- (diag(n)-phi[2]/lamb.max*C)
Q3 <- (diag(n)-phi[3]/lamb.max*C)

myrmvnorm <- function(n, mu, S) 
  sweep(matrix(rnorm(n*nrow(S)), n)%*%chol(S), 2, mu)

# defining the length of time series (number of years)
k <- 30

set.seed(1)

# simulating obsevational and innovation errors
w1 <- t(myrmvnorm(k, rep(0,n), solve(tau[1]*Q1)))
w2 <- t(myrmvnorm(k, rep(0,n), solve(tau[2]*Q2)))
w3 <- t(myrmvnorm(k, rep(0,n), solve(tau[3]*Q3)))

# generating the time series for observations and states
yy <- x1 <- x2 <- matrix(0, n, k)
x1[,1] <- w2[,1]
x2[,1] <- w3[,1]
for (i in 2:k) { 
  x2[,i] <- x2[,i-1] + w3[,i]
  x1[,i] <- x1[,i-1] + x2[,i-1] + w2[,i]
}
yy <- x1 + w1

# graphs of X1, X2 and w1 for area 10 and its neighbors
X11()
par(mfrow=c(3,1))
plot.ts(w1[10,])
for (i in 1:3)
  lines.ts(w1[nc.nb[[10]][i],], col=i)
plot.ts(x2[10,])
for (i in 1:3)
  lines.ts(x2[nc.nb[[10]][i],], col=i)
plot.ts(x1[10,])
for (i in 1:3)
  lines.ts(x1[nc.nb[[10]][i],], col=i, ylim=range(x1[nc.nb[[10]][i],]))


### defining the Cmatrix to use with model='generic1' for w1, w2 and w3
st.cmat <- kronecker(C, diag(k))
c.mat <- list(i=unlist(apply(st.cmat!=0, 1, which)),
              j=rep(1:nrow(st.cmat), rowSums(st.cmat!=0)), 
              values=st.cmat[st.cmat!=0])
##rm(st.cmat)


### building the augmented model
### ----------------------------
nd <- n*k
Y <- matrix(NA, nd*3-2*n, 3)
Y[1:nd             , 1] <- as.vector(t(yy))
Y[1:(nd-n) + nd    , 2] <- 0
Y[1:(nd-n) + 2*nd-n, 3] <- 0


### indices for the f() function
### ----------------------------
id1  <- (1:nd)[-((1:n)*k)]
id2  <- (1:nd)[-c(1,((1:(n-1))*k)+1)]
ix1  <- c(1:nd, id2, rep(NA,nd-n))                  ## indices for x1_t
ix1b <- c(rep(NA,nd), id1, rep(NA,nd-n))            ## indices for x1_{t-1}
wx1b <- c(rep(NA,nd), rep(-1,nd-n), rep(NA,nd-n))   ## weights for x1_{t-1}
ix2  <- c(rep(NA,nd),rep(NA,nd-n), id2)             ## indices for x2_t
ix2b <- c(rep(NA,nd), rep(id1, 2))                  ## indices for x2_{t-1}
wx2b <- c(rep(NA,nd), rep(-1,2*(nd-n)))             ## weights for x2_{t-1}
iw1  <- c(1:nd, rep(NA,2*(nd-n)))                   ## indices for w1_t
iw2  <- c(rep(NA,nd), id2, rep(NA,nd-n))            ## indices for w2_t
ww2  <- rep(c(NA,-1,NA), c(nd,length(id2),nd-n))    ## weights for w2_t
iw3  <- c(rep(NA,nd),rep(NA,nd-n), id2)             ## indices for w3_t
ww3  <- rep(c(NA,-1), c(nd+nd-n,length(id2)))       ## weights for w3_t


## formulating the model
## ---------------------
# with default priors for precision parameters and default initial values
formula6 = Y ~ f(iw1, model="generic1", Cmatrix=c.mat) +
               f(ix1, model="iid", initial=-10, fixed=T) +
               f(ix1b, wx1b, copy="ix1") +
               f(ix2, model="iid", initial=-10, fixed=T) +
               f(ix2b, wx2b, copy="ix2") +
               f(iw2, ww2, model="generic1", Cmatrix=c.mat) +
               f(iw3, ww3, model="generic1", Cmatrix=c.mat) -1

# with informative priors for precision parameters and default initial values
formula7 = Y ~ f(iw1, model="generic1", Cmatrix=c.mat, param=c(4,1/7.5, 0,0.001)) +
               f(ix1, model="iid", initial=-10, fixed=T) +
               f(ix1b, wx1b, copy="ix1") +
               f(ix2, model="iid", initial=-10, fixed=T) +
               f(ix2b, wx2b, copy="ix2") +
               f(iw2, ww2, model="generic1", Cmatrix=c.mat, param=c(4,0.08, 0,0.001)) +
               f(iw3, ww3, model="generic1", Cmatrix=c.mat, param=c(4,1/17.5, 0,0.001)) -1

# with vague priors for precision parameters and default initial values
formula8 = Y ~ f(iw1, model="generic1", Cmatrix=c.mat, param=c(0.01,1/3000, 0,0.001)) +
               f(ix1, model="iid", initial=-10, fixed=T) +
               f(ix1b, wx1b, copy="ix1") +
               f(ix2, model="iid", initial=-10, fixed=T) +
               f(ix2b, wx2b, copy="ix2") +
               f(iw2, ww2, model="generic1", Cmatrix=c.mat, param=c(0.01,1/5000, 0,0.001)) +
               f(iw3, ww3, model="generic1", Cmatrix=c.mat, param=c(0.01,1/7000, 0,0.001)) -1


# loading the INLA library
# ------------------------
require(INLA)


## call to fit the model
## ---------------------
r6 = inla(formula6, data = data.frame(ix1,ix1b,wx1b,ix2,ix2b,wx2b,iw1,iw2,iw3,ww2,ww3),
          family = rep("gaussian",3),
          control.family = list(list(initial=10, fixed=T),
          list(initial=10, fixed=T), list(initial=10, fixed=T)),
          control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

r7 = inla(formula7, data = data.frame(ix1,ix1b,wx1b,ix2,ix2b,wx2b,iw1,iw2,iw3,ww2,ww3),
          family = rep("gaussian",3),
          control.family = list(list(initial=10, fixed=T),
          list(initial=10, fixed=T), list(initial=10, fixed=T)),
          control.predictor=list(compute=TRUE, cdf=c(.025, .975)))

r8 = inla(formula8, data = data.frame(ix1,ix1b,wx1b,ix2,ix2b,wx2b,iw1,iw2,iw3,ww2,ww3),
          family = rep("gaussian",3),
          control.family = list(list(initial=10, fixed=T),
          list(initial=10, fixed=T), list(initial=10, fixed=T)),
          control.predictor=list(compute=TRUE, cdf=c(.025, .975)))


# elapsed  time
r6$cpu.used
r7$cpu.used
r8$cpu.used

par(mfrow=c(2,2))

# graph for observations (y)
rang <- range(r7[[17]][1:nd, 3:5], yy)
plot.ts(r7[[17]][1:nd,1], type="l", ylim=rang, col="red", ylab="y",xlab="time")
lines(r7[[17]][1:nd,3], col="blue", lty=3)
lines(r7[[17]][1:nd,5], col="blue", lty=3)
lines(as.vector(t(yy)))
legend("topleft",legend=c("observed","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("second order polynomial model")

#names(r7[[12]])

plot(as.vector(t(w1)), r7[[12]]$iw1$mean)
plot(as.vector(t(w2)), r7[[12]]$iw2$mean)
plot(as.vector(t(w3)), r7[[12]]$iw3$mean)
par(mfrow=c(2,2))
plot(as.vector(t(x1)), r7[[12]]$ix1$mean)
plot(as.vector(t(x2[,-1])), r7[[12]]$ix2$mean)


## graph of observations y_t for a given area and its neighbors
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
area = 20
for (i in c(area, nc.nb[[area]])) {
  id <- 1:k + (i-1)*k
  rang <- range(r7[[17]][id,3:5], yy[i,])
  plot.ts(r7[[17]][id,1], ylim=rang, col="red", ylab="y_t",xlab="time")
  lines(r7[[17]][id,3], col="blue", lty=3)
  lines(r7[[17]][id,5], col="blue", lty=3)
  lines(yy[i,])
  legend("topleft", legend=c("simulated y","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
}


## graph of X1_t state vector for a given area and its neighbors
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
area=20
for (i in c(area, nc.nb[[area]])) {
  id <- 1:k + (i-1)*k
  rang <- range(r7[[12]]$ix1[id, 4:6], x1[i,])
  plot.ts(r7[[12]]$ix1[id,2], ylim=rang, col="red", ylab="x1_t",xlab="time")
  lines(r7[[12]]$ix1[id,4], col="blue", lty=3)
  lines(r7[[12]]$ix1[id,6], col="blue", lty=3)
  lines(x1[i,])
  legend("topleft", legend=c("simulated X1","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
}

## graph of X2_t state vector for a given area and its neighbors
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,1,0))
area=20
for (i in c(area, nc.nb[[area]])) {
  id <- 1:(k-1) + (i-1)*(k-1)
  rang <- range(r7[[12]]$ix2[id, 4:6], x2[i,])
  plot.ts(r7[[12]]$ix2[id,2], ylim=rang, col="red", ylab="x2_t",xlab="time")
  lines(r7[[12]]$ix2[id,4], col="blue", lty=3)
  lines(r7[[12]]$ix2[id,6], col="blue", lty=3)
  lines(x2[i,2:k])
  legend("bottomleft", legend=c("simulated X2","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
}


# hyperparameter estimation values
# ---------------------------------
r6$summary.hyperpar
r7$summary.hyperpar
r8$summary.hyperpar


# posterior densities for precision parameters
# --------------------------------------------
par(mfrow=c(2,2))
rang1 <- range(r6[[23]][[1]][,2],r7[[23]][[1]][,2],r8[[23]][[1]][,2])
plot(r6[[23]][[1]], type="l", main=expression(tau[W1]), ylab="",xlab="" ,ylim=rang1)  # precision of V
  lines(r7[[23]][[1]],lty=3)
  lines(r8[[23]][[1]],lty=2)
abline(v=tau[1],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang3 <- range(r6[[23]][[3]][,2],r7[[23]][[3]][,2],r8[[23]][[3]][,2])
plot(r6[[23]][[3]], type="l", mai=expression(tau[W2]), ylab="" ,ylim=rang3,xlab="")  # precision of V
  lines(r7[[23]][[3]],lty=3)
  lines(r8[[23]][[3]],lty=2)
abline(v=tau[2],col=2, lwd=3)
legend("topright", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang5 <- range(r6[[23]][[5]][,2],r7[[23]][[5]][,2],r8[[23]][[5]][,2])
plot(r6[[23]][[5]], type="l", main=expression(tau[W3]), ylab="" ,ylim=rang5,xlab="")  # precision of V
  lines(r7[[23]][[5]],lty=3)
  lines(r8[[23]][[5]],lty=2)
abline(v=tau[3],col=2, lwd=3)
legend("topright", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)


par(mfrow=c(2,2))
rang2 <- range(r6[[23]][[2]][,2],r7[[23]][[2]][,2],r8[[23]][[2]][,2])
plot(r6[[23]][[2]], type="l", main=expression(phi[W1]), ylab="" ,ylim=rang2,xlab="")  
  lines(r7[[23]][[2]],lty=3)
  lines(r8[[23]][[2]],lty=2)
abline(v=phi[1],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang4 <- range(r6[[23]][[4]][,2],r7[[23]][[4]][,2],r8[[23]][[4]][,2])
plot(r6[[23]][[4]], type="l", main=expression(phi[W2]), ylab="" ,ylim=rang4,xlab="")  
  lines(r7[[23]][[4]],lty=3)
  lines(r8[[23]][[4]],lty=2)
abline(v=phi[2],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)
rang6 <- range(r6[[23]][[6]][,2],r7[[23]][[6]][,2],r8[[23]][[6]][,2])
plot(r6[[23]][[6]], type="l", main=expression(phi[W3]), ylab="" ,ylim=rang6,xlab="")  
  lines(r7[[23]][[6]],lty=3)
  lines(r8[[23]][[6]],lty=2)
abline(v=phi[3],col=2, lwd=3)
legend("topleft", legend=c("with vague prior","with informative prior","with default prior"), col=c("black", "black", "black"),lty=c(3,1,2),bty="n",cex=1)


require(RColorBrewer)
c4 <- c(rev(brewer.pal(3, "Blues"))[1:2], brewer.pal(3, "Reds")[2:3])

### map of y at times 2, 7 and 15
postscript("D:\\meus artigos\\INLA-DLMs\\figures/ex6_mapy_2715-c.eps", 
  horizontal=FALSE, width=17, height=5)
par(mfrow=c(3,2), mar=c(0,0,0,0))
ind1 = c(2,7,15)
for (t in ind1) {
  q <- quantile(c(yy[,t], r7[[17]][seq(1, nd, k)+t-1,1]), 0:4/4)
  plot(nc, col=c4[findInterval(yy[,t], q, T)])
#title(t, "areas")
  plot(nc, col=c4[findInterval(r7[[17]][seq(1, nd, k)+t-1,1], q, T)])
#title("areas")
  legend("bottomleft", leglabs(format(q,dig=1),bet="to"), fill=c4, bty="n", ncol=2)
}
dev.off()

### map of x1 at times 2, 7 and 15
postscript("D:\\meus artigos\\INLA-DLMs\\figures/ex6_mapx1_2715-c.eps", 
  horizontal=FALSE, width=17, height=5)
par(mfrow=c(3,2), mar=c(0,0,0,0))
ind1 = c(2,7,15)
for (t in ind1) {
  q <- quantile(c(x1[,t], r7[[12]]$ix1$mean[seq(1,nd,k)+t-1]), 0:4/4)
  plot(nc, col=c4[findInterval(x1[,t], q, T)])
  plot(nc, col=c4[findInterval(r7[[12]]$ix1$mean[seq(1, nd, k)+t-1], q, T)])
  legend("bottomleft", leglabs(format(q,dig=1),bet="to"), fill=c4, bty="n", ncol=2)
}
dev.off()

### map of x2 at times 2, 7 and 15
postscript("D:\\meus artigos\\INLA-DLMs\\figures/ex6_mapx2_2715-c.eps", 
  horizontal=FALSE, width=17, height=5)
par(mfrow=c(3,2), mar=c(0,0,0,0))
ind1 = c(2,7,15)
for (t in ind1) {
  q <- quantile(c(x2[,t], r7[[12]]$ix2$mean[seq(1,nd-n,k-1)+t-2]), 0:4/4)
  plot(nc, col=c4[findInterval(x2[,t], q, T)])
  plot(nc, col=c4[findInterval(r7[[12]]$ix2$mean[seq(1, nd-n, k-1)+t-2], q, T)])
  legend("bottomleft", leglabs(format(q,dig=1),bet="to"), fill=c4, bty="n", ncol=2)
}
dev.off()



########################
#    CASE STUDIES      #
########################

## -------------------------------------------------------------------------------------------------
## Example 7: quarterly UK gas consumption from 1960 to 1986, in millions of therms 
## (for details see Durbin and Koopman 2001, p. 233).
## model with a first order polynomial trend with time-varying coefficients and an
## unstructured seasonal component, also varying over time.
## The response is the (base 10) logarithm of the UK gas consumption (assumed as normal distributed)
##
## Observational equation:
##    y_t = T_t + S_t + v_t  ,   t=1,...,n                   (1)
##
## System equations:
##    T_t     = T_{t-1} + \beta_{t-1} + w_{1,t} ,            (2)    t=2,...,n
##    \beta_t = \beta_{t-1} + w_{2,t} ,                      (3)    t=2,...,n
##    S_t     = -(S_{t-1} + S_{t-2} + S_{t-3}) + w_{3,t} ,   (4)    t=4,...,n
##
## We can re-write Eq. (2) as
##    0 = T_t - T_{t-1} - \beta_{t-1} + w_{1,t} ,                   t=2,...,n
##
## and then merge it with observational equation (1) in an augmented model with two
## different likelihoods. The first n datapoints being Gaussian distributed with
## unknown precision, whereas the last `n-1' datapoints, which are forced
## to be 0, are observed with a high and fixed precision. 
##
## For the slope term in Eq. (3) we use an "RW1" model.
## For the seasonal term in Eq. (4) we use a "seasonal" model.
## ---------------------------------------------------------------------------------------------------

# loading the data set
# --------------------
data(UKgas)
plot(log(UKgas))

n = length(UKgas)
m = n-1


# building the augmented model
# ----------------------------
Y = matrix(NA, n+m, 2)
Y[1:n,             1] = log10(UKgas)
Y[1:m + n,         2] = 0


## indices for the INLA library
# ----------------------------
i       = c(1:n, 2:n)                # T_t
j       = c(rep(NA,n), 1:m)          # T_{t-1}
weight1 = c(rep(NA,n), rep(-1,m))    # weights for j
l       = c(rep(NA,n), 1:m)          # \beta_{t-1}
weight2 = c(rep(NA,n), rep(-1,m))    # weights for l
w1      = c(rep(NA,n), 2:n)          # w_{1,t}
q       = c(1:n, rep(NA,m))          # S_t 


# formulating the model
# ---------------------
formula = Y ~ f(q, model="seasonal", season.length=4,param=c(1,0.0001),initial=4) +
              f(l, weight2, model="rw1",param=c(1,0.0001),initial=4, constr=F) +
              f(i, model="iid", initial=-10, fixed=TRUE) +
              f(j, weight1, copy="i") + 
              f(w1, model ="iid") -1


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r = inla(formula, data = data.frame(i,j,weight1,l,weight2,q,w1),
         family = rep("gaussian",2),
         control.family = list(list(),list(initial=10, fixed=TRUE)),
         control.predictor=list(compute=TRUE))


# elapsed time
r$cpu.used


## plotting the results
# ---------------------
rang <- range(r[[17]][1:n, 3:5], log10(UKgas))
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="z",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(log10(UKgas)[1:n])
legend("topleft",legend=c("observed","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("Trend+Seasonal model")


# fitting the data using the sspir package (Dethlefsen and Lundbye-Christensen, 2006)
# (it uses an iterated extended Kalman filtering approach)
# --------------------------------------------------------
require(sspir)
phistart <- StructTS(log10(UKgas),type="BSM")$coef[c(4,1,2,3)]
gasmodel <- ssm(log10(UKgas) ~ -1 + tvar(polytime(time,1)) +
                tvar(sumseason(time,4)), phi=phistart)
fit <- getFit(gasmodel)


# Decomposition of the UK gas time series with INLA and sspir
# based on code by Dethlefsen and Lundbye-Christensen (JSS, 16, 2006)
# -------------------------------------------------------------------
par(mfrow=c(3,1),mar=c(0,5.1,0,2.1),oma=c(6,0,5,0))
plot.ts(r[[12]]$i$mean,ylab="Trend",axes=FALSE,col=4)
lines(r[[12]]$i[,4],lty=3,col=4)
lines(r[[12]]$i[,6],lty=3,col=4)
lines(as.vector(fit$m[,1]), col=2)
#abline(v=44,lty=2)
#abline(v=80,lty=2)
axis(2)
box()

rang <- range(r[[12]]$l$mean,r[[12]]$l[,4],r[[12]]$l[,6])
plot.ts(r[[12]]$l$mean,ylab="Slope",axes=FALSE,col=4,lwd=2,ylim=rang)
lines(r[[12]]$l[,4],col=4,lty=3)
lines(r[[12]]$l[,6],col=4,lty=3)
lines(as.vector(fit$m[1:107,2]), col=2, lwd=2)
#abline(v=43,lty=2)
#abline(v=79,lty=2)
axis(2)
box()

plot.ts(r[[12]]$q$mean,ylab="Season",xlab="Time",axes=FALSE,col=4)
lines(r[[12]]$q[,4],lty=3,col=4)
lines(r[[12]]$q[,6],lty=3,col=4)
lines(as.vector(fit$m[,3]), col=2)
abline(v=44,lty=2)
abline(v=80,lty=2)
axis(2)
axis(1, seq(0,108,4), 1960:1987, las=3)
box()


# hyperparameter estimation values
# --------------------------------
r$summary.hyperpar


# posterior density for precision parameters
# ------------------------------------------
par(mfrow=c(2,2))
plot(r[[23]][[1]], type="l", xlab="1/V", main="precision for observations", ylab="")  # precision of V
plot(r[[23]][[2]], type="l", xlab="1/W1", main="precision for w_1t", ylab="")  # precision of W
plot(r[[23]][[3]], type="l", xlab="1/W2", main="precision for w_2t", ylab="")  # precision of W
plot(r[[23]][[4]], type="l", xlab="1/W3", main="precision for w_3t", ylab="")  # precision of W




## ----------------------------------------------------------------------------------------------
## Example 8: monthly numbers of light goods van drivers killed in road accidents, from
## January 1969 to December 1984 (for details see Harvey and Durbin, 1986).
## This is an state space model for Poisson data with a 13-dimensional latent process, consisting
## of an intervention parameter, seatbelt, changing value from zero to one in February
## 1983, a constant monthly seasonal, and a trend modelled as a random walk.
##
## Observational equation:
##    y_t \sim Poisson(\mu_t)
##    log(\mu_t) = \lambda_t = T_t + \alpha*seatbelt + S_t  ,   t=1,...,n
##
## System equations:
##    T_t     = T_{t-1} + w_{1,t}      ,                        t=2,...,n
##    S_t     = -(S_{t-1} + ... + S_{t-11})                     t=12,...,n
##
## -----------------------------------------------------------------------------------------------


# loading the data set
# --------------------
require(sspir)
data("vandrivers")
plot(vandrivers$y, type="l")

set.seed(123456)

n = dim(vandrivers)[1]
y = vandrivers$y
belt = vandrivers$seatbelt
i = j = 1:n                 # indices for T_t and S_t


# formulating the model
# ---------------------
formula = y ~ belt + f(i, model="rw1",param=c(1,0.0005), constr=F) +  
              f(j, model="seasonal", season.length=12) -1


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r = inla(formula, data = data.frame(belt,i,j,y),
         family = "poisson",
         control.predictor=list(compute=TRUE))


## plotting the results
rang <- range(r[[17]][1:n, 3:5], vandrivers$y)
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="vandrivers",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(vandrivers$y[1:n])
legend("topright", legend=c("observed","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")
#title("Van drivers example")


# fitting the data using the sspir package (Dethlefsen and Lundbye-Christensen, 2006)
# (it uses an iterated extended Kalman filtering approach)
# --------------------------------------------------------
vd <- ssm( y ~ tvar(1) + seatbelt + sumseason(time,12),
          family=poisson(link="log"), data=vandrivers,
          phi=c(1,0.0004945505),
          C0=diag(13)*1000)
vd.fit <- getFit(vd)


## plotting the trend + intervention effect with the two approaches
# -----------------------------------------------------------------
plot(1:192, vandrivers$y,ylim=range(vandrivers$y),type="p",ylab="vandrivers killed",xlab="time",axes=FALSE)
lines(1:192, exp(vandrivers$seatbelt*r[[2]][1] + r[[12]][[1]][,2]),col="blue")
lines(1:192, exp(vandrivers$seatbelt*r[[2]][1] + r[[12]][[1]][,4]),lty=3,col="blue")
lines(1:192, exp(vandrivers$seatbelt*r[[2]][1] + r[[12]][[1]][,6]),lty=3,col="blue")
lines(1:192, exp(vandrivers$seatbelt*vd.fit$m[,2] + vd.fit$m[,1]),col=2 )
axis(2)
axis(1, seq(0,192,12), 1969:1985, las=3)
box()


# effect of seat belt legislation (% of reduction in deaths)
100 * (1 - exp(r[[2]][1]))           # with INLA
100 * (1 - exp(vd.fit$m[1,2]))       # with sspir package


# summary of fixed effects
# ------------------------
r$summary.fixed


# summary of posterior hyperparameters
# ------------------------------------
r$summary.hyperpar


# posterior density for precision parameter
# -----------------------------------------
par(mfrow=c(2,1))
plot(r[[23]][[1]], type="l", xlab="1/V", main="precision for w_t", ylab="")  # precision for w_t


# posterior density for the intervention effect
# ---------------------------------------------
plot(r[[3]][[1]], type="l", xlab="intervention effect", main="", ylab="")  # 
  abline(v=r[[2]][1,1],col=2, lwd=3)
  abline(v=r[[2]][1,3],col=4, lty=2)
  abline(v=r[[2]][1,5],col=4, lty=2)




## ----------------------------------------------------------------------------------------------
## Example 9: Monthly registered cases of mumps in New York City, January 1928 through June 1972 
## (for details see Hipel and McLeod, 1994 and Dethlefsen and Lundbye-Christensen, 2006). 
## Mumps incidence is modelled with a first order polynomial trend with time-varying coefficients
## and a time-varying harmonic seasonal component.
##
## Observational equation:
##    y_t \sim Poisson(\mu_t)
##    log(\mu_t) = \lambda_t = T_t + H_t  ,   t=1,...,n                     (1)
##
## System equations:
##    T_t     = T_{t-1} + \beta_{t-1} + w_{1,t}  ,              t=2,...,n   (2)
##    \beta_t = \beta_{t-1} + w_{2,t}                           t=2,...,n   (3)
##    H_t     = a_t*cos((2\pi/12)*t) + b_t*sin((2\pi/12)*t)     t=1,...,n   
##    a_t = a_{t-1}*cos((2\pi/12)*t) + w_{3,t}                  t=1,...,n   (4)
##    b_t = b_{t-1}*sin((2\pi/12)*t) + w_{4,t}                  t=1,...,n   (5)
##
## We can re-write Eq. (2) as
##    0 = T_t - T_{t-1} - \beta_{t-1} + w_{1,t}  ,              t=2,...,n
##
## and then merge it with observational equation (1) in an augmented model with two
## different likelihoods. The first n datapoints being Poisson distributed with
## unknown precision, whereas the last `n-1' datapoints, which are forced
## to be 0, are observed with a high and fixed precision. 
##
## For \beta_t and the seasonal terms a_t and b_t in Eq. (3) to (5) we use an RW1 model.
## ----------------------------------------------------------------------------------------------

# loading the data set
# --------------------
require(sspir)
data("mumps")
plot(mumps,type="l")

n = length(mumps)
m = n-1
t = 1:n
cosw = cos((2*pi/12)*t)
sinw = sin((2*pi/12)*t)

set.seed(123456)

# building the augmented model
# ----------------------------
Y = matrix(NA, n+m, 2)
Y[1:n,             1] = mumps
Y[1:m + n,         2] = 0


## indices for the INLA library
# ----------------------------
i       = c(1:n, 2:n)                # T_t
j       = c(rep(NA,n), 1:m)          # T_{t-1}
weight1 = c(rep(NA,n), rep(-1,m))    # weights for j
l       = c(rep(NA,n), 1:m)          # \beta_{t-1}
weight2 = c(rep(NA,n), rep(-1,m))    # weights for l
w1      = c(rep(NA,n), 2:n)          # w_{1,t}
q       = c(1:n, rep(NA,m))          # a_t 
cosine  = c(cosw,rep(NA,m))          # weights for a_t
rr      = c(1:n, rep(NA,m))          # b_t 
sine    = c(sinw,rep(NA,m))          # weights for b_t


# formulating the model
# ---------------------
formula = Y ~ f(q, cosine, model="rw1",param=c(1,0.01),initial=4, constr=F) +
              f(rr, sine, model="rw1",param=c(1,0.01),initial=4, constr=F) + 
              f(l, weight2, model="rw1",param=c(1,0.2),initial=4, constr=F) +
              f(i, model="iid", initial=-10, fixed=TRUE) +
              f(j, weight1, copy="i") + 
              f(w1, model ="iid") -1


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r = inla(formula, data = data.frame(cosine,sine,i,j,weight1,l,weight2,q,rr,w1),
         family = c("poisson","gaussian"),
         control.family = list(list(),list(initial=10, fixed=TRUE)),
         control.predictor=list(compute=TRUE))


# fitting the data using the sspir package (Dethlefsen and Lundbye-Christensen, 2006)
# (it uses an iterated extended Kalman filtering approach)
# --------------------------------------------------------
index <- 1:length(mumps)
mumps.m <- ssm( mumps ~ -1 + tvar(polytime(index,1)) +
                tvar(polytrig(index,12,1)), family=poisson(link=log),
                phi=c(0,0,0.0005,0.0001),
                C0 = diag(4))
m3.fit <- getFit(mumps.m)


# comparison with sspir package
# -----------------------------
par(mfrow=c(2,2))

plot(1:534, m3.fit$m[,1], type="l")
lines(r[[12]]$i$mean, col=2)

plot(1:534, m3.fit$m[,2], type="l")
lines(r[[12]]$l$mean, col=2)

plot(1:534, m3.fit$m[,3], type="l")
lines(r[[12]]$q$mean, col=2)

plot(1:534, m3.fit$m[,4], type="l")
lines(r[[12]]$rr$mean, col=2)


# Plotting the results of variation in the incidence of mumps cases
# based on code by Dethlefsen and Lundbye-Christensen (JSS, 16, 2006)
# -------------------------------------------------------------------
par(mfrow=c(3,1),mar=c(0,5.1,0,2.1),oma=c(6,0,5,0))
plot(1:534,mumps,type="l",axes=FALSE)
lines(1:534,exp(r[[12]]$i$mean),type="l",col=4,lwd=2)
lines(1:534,exp(m3.fit$m[,1]),type='l',col=2,lwd=2)
axis(2)
box()

plot(1:534,12*atan2(m3.fit$m[,4],m3.fit$m[,3])/(2*pi),type='l',ylim=c(2.7,5.3),ylab='Peak',xlab='',lwd=2,axes=FALSE,col=2)
lines(1:534,12*atan2(r[[12]]$rr$mean,r[[12]]$q$mean)/(2*pi),type='l',ylim=c(2.7,5.3),ylab='Peak',xlab='',lwd=2,axes=FALSE,col=4)
abline(h=3,lty=3)
abline(h=4,lty=3)
abline(h=5,lty=3)
axis(2,at=c(3,4,5),labels=c("Apr","May","Jun"))
box()

plot(1:534,exp(2*sqrt(m3.fit$m[,3]^2 + m3.fit$m[,4]^2)),type='l',ylab='PT-ratio',xlab='Years',ylim=c(0,12),lwd=2,axes=FALSE,col=2)
lines(1:534,exp(2*sqrt((r[[12]]$q$mean)^2 + (r[[12]]$rr$mean)^2)),type='l',ylab='PT-ratio',xlab='Years',ylim=c(0,12),lwd=2,axes=FALSE,col=4)
abline(h=0,lty=3)
abline(h=5,lty=3)
abline(h=10,lty=3)
axis(2)
axis(1, seq(0,534,12), 1928:1972, las=3)
box()


## observed vs predicted values
# -----------------------------
rang <- range(r[[17]][1:n, 3:5], mumps)
plot(r[[17]][1:n,1], type="l", 
 ylim=rang, col="red", xlim=c(1,n),ylab="mumps cases",xlab="time")
lines(r[[17]][1:n,3], col="blue", lty=3)
lines(r[[17]][1:n,5], col="blue", lty=3)
lines(mumps)
legend("topright", legend=c("observed","posterior mean","95% CI"), col=c("black", "red","blue"),lty=c(1,1,2),bty="n")




## -----------------------------------------------------------------------------------------------------------
## Example 10: A dynamic regression model with three covariates to analyse market share 
## for a consumer product. Data are registered weekly for 1990 and 1991 
## (for details see Pole, West and Harrison, 1994, ch5). 
##
## Observational equation:
##    y_t = \beta_{0,t} + \beta_{1,t}*PRICE_t + \beta_{2,t}*PROM_t + \beta_{3,t}*CPROM_t + v_t  ,   t=1,...,n
##
## System equations:
##    \beta_{1,t} = \beta_{1,t-1} + w_{2,t}  ,  t=2,...,n
##    \beta_{2,t} = \beta_{2,t-1} + w_{3,t}  ,  t=2,...,n
##    \beta_{3,t} = \beta_{3,t-1} + w_{4,t}  ,  t=2,...,n
##
## and \beta_0 is a fixed intercept
## -----------------------------------------------------------------------------------------------------------

# loading the data set
share0 <- read.table("share.txt", header=TRUE)
share <- share0
share[34,4] <- NA     # removing an outlier observation at week 34/1990

# loading BATS results
bats=read.table("BATS_share_output.txt",header=T)


# defining indices for rw1 coefficients
# -------------------------------------
n <- nrow(share)
id3 <- id2 <- id1 <- id <- 1:n


# formulating the model
# ---------------------
formula <- SHARE ~ f(id1, PROM, model="rw1", param=c(1,0.001), constr=F)  +
                   f(id2, PRICE, model="rw1", param=c(1,0.01), constr=F) + 
                   f(id3, CPROM, model="rw1", param=c(1,0.001), constr=F) 


# loading the INLA library
# ------------------------
require(INLA)


# call to fit the model
# ---------------------
r.rw1 <- inla(formula, data=data.frame(id,id1,id2,id3,share),
              quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
              control.predictor=list(compute=TRUE))


# elapsed time
r.rw1$cpu.used

par(mfrow=c(2,2))
for (i in 1:4) plot(r.rw1$marginals.hyperpar[[i]], type="l")

# summary for hyperparameters
r.rw1$summary.hyperpar


## plotting the results
# ---------------------

# regression coefficients
par(mfrow=c(3,1),mar=c(0,5.1,0,2.1),oma=c(6,0,5,0))
ylm <- range(r.rw1[[12]][[1]][,4:8]) #,bats$PROM)
plot.ts(r.rw1[[12]][[1]]$mean,type='l',
        ylab='PROM',xlab='',axes=FALSE, col=2,ylim=ylm)
lines.ts(r.rw1[[12]][[1]][,4], col=4, lty=3)
lines.ts(r.rw1[[12]][[1]][,8], col=4, lty=3)
#lines.ts(bats$PROM)
axis(2)
box()

ylm <- range(r.rw1[[12]][[2]][,4:8]) #,bats$PRICE)
plot.ts(r.rw1[[12]][[2]]$mean,type='l',
        ylab='PRICE',xlab='',axes=FALSE, col=2,ylim=ylm)
lines.ts(r.rw1[[12]][[2]][,4], col=4, lty=3)
lines.ts(r.rw1[[12]][[2]][,8], col=4, lty=3)
#lines.ts(bats$PRICE)
axis(2)
box()

ylm <- range(r.rw1[[12]][[3]][,4:8]) #,bats$CPROM)
plot.ts(r.rw1[[12]][[3]]$mean,type='l',
        ylab='CPROM',axes=FALSE, col=2,ylim=ylm,xlab="Time")
lines.ts(r.rw1[[12]][[3]][,4], col=4, lty=3)
lines.ts(r.rw1[[12]][[3]][,8], col=4, lty=3)
#lines.ts(bats$CPROM)
axis(2)
axis(1, seq(0,104,52), 1990:1992, las=3)
#axis(1)
box()



# predicted values and estimated level with INLA
par(mfrow=c(1,1))
with(share0, plot.ts(SHARE, type="p",ylim=c(38,45),ylab="% Market share",axes=FALSE))
#lines(fsh$mu, col=4)
lines(r.rw1[[17]]$mean, col=2)
lines(r.rw1[[17]][,3], col=4, lty=3)
lines(r.rw1[[17]][,7], col=4, lty=3)
abline(h=r.rw1$summary.fixed[1,1])
abline(h=r.rw1$summary.fixed[1,4],lty=3)
abline(h=r.rw1$summary.fixed[1,6],lty=3)
legend("bottomright", legend=c("posterior mean","90% CI"), col=c("red","blue"),lty=c(1,2),bty="n")
axis(2)
axis(1, seq(0,104,52), 1990:1992, las=3)
box()

# predicted values and estimated level with BATS
par(mfrow=c(1,1))
with(share0, plot.ts(SHARE, type="p",ylim=c(38,45),ylab="% Market share",axes=FALSE))
#lines(fsh$mu, col=4)
lines(bats$forecast, col=2)
lines(bats$forecast-bats$X90CI, col=4, lty=3)
lines(bats$forecast+bats$X90CI, col=4, lty=3)
abline(h=41.68)
abline(h=41.51,lty=3)
abline(h=41.85,lty=3)
legend("bottomright", legend=c("posterior mean","90% CI"), col=c("red","blue"),lty=c(1,2),bty="n")
axis(2)
axis(1, seq(0,104,52), 1990:1992, las=3)
box()


# Forecasts (five steps head)
# ---------------------------

# defining indices for rw1 coefficients
# -------------------------------------
m <- n + 5
id3f <- id2f <- id1f <- idf <- 1:m

SHAREf <- c(share$SHARE, rep(NA,5))

# scenario 1:
PROMf = c(share$PROM, rep(0,5))
CPROMf = c(share$CPROM, rep(0,5))
PRICEf = c(share$PRICE, rep(share$PRICE[n],5))

# formulating the model
# ---------------------
formula <- SHAREf ~ f(id1f, PROMf, model="rw1", param=c(2,0.01), constr=F)  +
                    f(id2f, PRICEf, model="rw1", param=c(2,0.01), constr=F) + 
                    f(id3f, CPROMf, model="rw1", param=c(2,0.01), constr=F)  

# call to fit the model
# ---------------------
r.rw1.f1 <- inla(formula, data=data.frame(idf,id1f,id2f,id3f,PROMf,CPROMf,PRICEf,SHAREf),
                 quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 control.predictor=list(compute=TRUE))

r.rw1.f1[[17]][[1]][105:109]
r.rw1.f1[[17]][[4]][105:109]
r.rw1.f1[[17]][[6]][105:109]

# scenario 2:
PROMf = c(share$PROM, share$PROM[1:5])
CPROMf = c(share$CPROM, rep(0,5))
PRICEf = c(share$PRICE, rep(share$PRICE[n],5))

# formulating the model
# ---------------------
formula <- SHAREf ~ f(id1f, PROMf, model="rw1", param=c(1,0.001), constr=F)  +
                    f(id2f, PRICEf, model="rw1", param=c(1,0.01), constr=F) + 
                    f(id3f, CPROMf, model="rw1", param=c(1,0.001), constr=F) 

# call to fit the model
# ---------------------
r.rw1.f2 <- inla(formula, data=data.frame(idf,id1f,id2f,id3f,PROMf,CPROMf,PRICEf,SHAREf),
                 quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 control.predictor=list(compute=TRUE))

r.rw1.f2[[17]][[1]][105:109]
r.rw1.f2[[17]][[4]][105:109]
r.rw1.f2[[17]][[6]][105:109]

# scenario 3:
PROMf = c(share$PROM, rep(0,5))
CPROMf = c(share$CPROM, share$CPROM[1:5])
PRICEf = c(share$PRICE, rep(share$PRICE[n],5))

# formulating the model
# ---------------------
formula <- SHAREf ~ f(id1f, PROMf, model="rw1", param=c(1,0.001), constr=F)  +
                    f(id2f, PRICEf, model="rw1", param=c(1,0.01), constr=F) + 
                    f(id3f, CPROMf, model="rw1", param=c(1,0.001), constr=F) 

# call to fit the model
# ---------------------
r.rw1.f3 <- inla(formula, data=data.frame(idf,id1f,id2f,id3f,PROMf,CPROMf,PRICEf,SHAREf),
                 quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 control.predictor=list(compute=TRUE))

r.rw1.f3[[17]][[1]][105:109]
r.rw1.f3[[17]][[4]][105:109]
r.rw1.f3[[17]][[6]][105:109]

# scenario 4:
PROMf = c(share$PROM, share$PROM[1:5])
CPROMf = c(share$CPROM, share$CPROM[1:5])
PRICEf = c(share$PRICE, rep(share$PRICE[n],5))

# formulating the model
# ---------------------
formula <- SHAREf ~ f(id1f, PROMf, model="rw1", param=c(1,0.001), constr=F)  +
                    f(id2f, PRICEf, model="rw1", param=c(1,0.01), constr=F) + 
                    f(id3f, CPROMf, model="rw1", param=c(1,0.001), constr=F) 

# call to fit the model
# ---------------------
r.rw1.f4 <- inla(formula, data=data.frame(idf,id1f,id2f,id3f,PROMf,CPROMf,PRICEf,SHAREf),
                 quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 control.predictor=list(compute=TRUE))

r.rw1.f4[[17]][[1]][105:109]
r.rw1.f4[[17]][[4]][105:109]
r.rw1.f4[[17]][[6]][105:109]

