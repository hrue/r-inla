gen.cure.gompertz = function(n,a,b,p) {
    if (length(b) == 1) b <- rep(b, n)
    if (length(p) == 1) p <- rep(p, n)
    rm = rbinom(n=n,size=1,prob=1-p)
    t=rep(NA,n)
    for(i in 1:n){
        t[i]=ifelse(rm[i]==0,
                    Inf, 
                    log((-(a/b[i])*log(1-runif(n=1,min=0,max=1-p[i]))) + 1)*(1/a))
    }
    t_finite = ifelse(t==Inf,0,t)
    u2 = runif(n=n,0,max(t_finite))
    t2 = pmin(t,u2)
    delta = ifelse(t<u2,1,0)
    return(cbind(t2,delta))
} 

n = 1000
a_par=-0.5 ##defective
x <- rnorm(n)
b_par <- exp(0.2 + x)
p_par = exp(b_par/a_par) ##cure fraction

data = gen.cure.gompertz(n=n,a=a_par,b=b_par,p=p_par)
colnames(data) = c("tempo", "delta")
data <- as.data.frame(data)

data$x <- x
r <- inla(inla.surv(tempo, delta) ~ 1 + x,
          data = data, 
          family = "dgompertzsurv",
          control.inla = list(cmin = 0), 
          verbose = TRUE)
summary(r)
