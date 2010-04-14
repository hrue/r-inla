data(Kidney)
n = length(Kidney$event)

## make the input to the survival-module in inla()
surv.time = list(truncation=rep(0,n), event = Kidney$event,
                 lower=Kidney$time, upper = rep(0,n), time=Kidney$time)

formula =  surv.time ~ age+sex+dis1+dis2+dis3+f(ID,param=c(1,0.01),initial=0, model="iid")

d = c(as.list(Kidney), surv.time = surv.time)

model1=inla(formula,family="weibull",data=d, control.fixed=list(prec=1e-5),
            control.data=list(param=c(1,1)), keep = TRUE, verbose=TRUE)

summary(model1)
plot(model1)


