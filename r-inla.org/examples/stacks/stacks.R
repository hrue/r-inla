# read the data and set names
data = as.data.frame(read.table("stacks.dat"))
names(data) = c("y","air","temp","acid")

# set priors
param.data = list(prec = list(param = c(1.0e-3, 1.0e-3)))

# this is the simple regression model with gaussian data

formula1 = y ~  air + temp + acid
mod1 = inla(formula1, data=data,
        control.fixed = list(prec = 0.00001),
        control.family = list(hyper = param.data))


# for the ridge regression the three parameters beta1..beta3 have a
# common unknown variance which we have to estimate

# to implement this we have to use copies, and change the data set
data$beta1 = rep(1,21)
data$beta2 = rep(2,21)
data$beta3 = rep(3,21)

# this is the prior for the precision of beta
param.beta = list(prec = list(param = c(1.0e-3, 1.0e-3)))

formula.ridge  = y  ~ f(beta1, air, model="iid", values = c(1,2,3),
        hyper = param.beta) +
    f(beta2, temp, copy="beta1", fixed=T) +
    f(beta3, acid, copy="beta1", fixed=T)


mod.ridge = inla(formula.ridge, data=data)


# compare results
simple = mod1$summary.fixed[,c(1,2)]
ridge = rbind(mod.ridge$summary.fixed[,c(1,2)], mod.ridge$summary.random$beta1[,c(2,3)])
colnames(ridge) = colnames(simple)
rownames(ridge) = rownames(simple)

print(simple)
print(ridge)



