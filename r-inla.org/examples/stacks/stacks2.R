library(INLA)
# read the data and set names
data = as.data.frame(read.table("stacks.dat"))
names(data) = c("y","air","temp","acid")
n = nrow(data)

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
data$beta1 = rep(1,n)
data$beta2 = rep(2,n)
data$beta3 = rep(3,n)

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

# A-matrix method: rigde regression as a mixed model
# y = X beta + Z b + e
# (beta = fixed effects == intercept,  b = random effects for coefficients of predictors)
X = matrix(1,nrow = n, ncol= 1)           # intercept
Z = as.matrix(data[,1 + 1:3]) # predictors

pX = ncol(X); pZ = ncol(Z)

idx.X = c(1:pX, rep(NA, pZ))
idx.Z = c(rep(NA,pX), 1:pZ)

hyper.fixed = list(prec = list(initial = log(0.00001), fixed=TRUE))
param.Z =  param.beta #  param.beta = list(prec = list(param = c(1.0e-3, 1.0e-3))) # defined earlier
formula = y ~ -1 + f(idx.X,  model="iid", hyper = hyper.fixed) + f(idx.Z,  model="iid", hyper = param.Z)
result1 = inla(formula, data = list(y=data$y, idx.X=idx.X, idx.Z=idx.Z),
               control.predictor = list(A=cbind(X, Z)), control.family = list(hyper = param.data))   
# explicit precision matrix for coefficients of Z, in this case: identity
Qz = diag(pZ)
formula = y ~ -1 + f(idx.X,  model="iid", hyper = hyper.fixed) + f(idx.Z,  model="generic", Cmatrix=Qz, hyper = param.Z)
result2 = inla(formula, data = list(y=data$y, idx.X=idx.X, idx.Z=idx.Z),
              control.predictor = list(A=cbind(X, Z)), control.family = list(hyper = param.data))

all.equal(result1$summary.random$idx.Z[,-1] , result2$summary.random$idx.Z[,-1]) # TRUE

ridge2 = rbind(result1$summary.random$idx.X[,c(2,3)], result1$summary.random$idx.Z[,c(2,3)])
colnames(ridge2) = colnames(simple)
rownames(ridge2) = rownames(simple)
print(ridge)
print(ridge2)






