remove(list=ls())
require(INLA)
inla.my.update()

n = 50000
factor = 0.0001

field = inla.spde.make.index("field", n.mesh=n)
covar1 = matrix(runif(n),n,1)
covar2 = matrix(runif(n),n,1)
A = inla.spde.make.A(n.mesh=n)
y = covar1+covar2+rnorm(n)*1

formula = y ~ -1 + covar1+covar2 + f(field, model="rw1")

stack = list()
stack[[1]] =
    inla.stack(data = list(y=y),
               A = list(A, 1),
               effect =
               list(field,
                    list(covar1=covar1, covar2=covar2)))
stack[[2]] =
    inla.stack(data = list(y=y),
               A = list(A, covar1, covar2),
               effect =
               list(field,
                    covar1=1,
                    covar2=1))

if (FALSE) {
A.covar = cBind(covar1, covar2)
Q0 = Diagonal(n)
Qc = Diagonal(2)
cBind(Q0+t(A.covar)%*%Qc%*%A.covar, -A.covar)
cBind(-t(A.covar), Qc)
}

result = list()
result. = list()
for (k in 1:length(stack)) {
##    result[[k]] =
##        inla(
##             formula, data=inla.stack.data(stack[[k]]), family="gaussian",
##             control.predictor=list(A=inla.stack.A(stack[[k]])),
##             verbose=TRUE
##         )
    result.[[k]] =
        inla(
             formula, data=inla.stack.data(stack[[k]]), family="gaussian",
             control.predictor=list(A=inla.stack.A(stack[[k]])),
             control.compute=list(q=TRUE),
             control.inla=list(
##             global.node.factor=factor
             reordering="AMD"
             ),
             verbose=TRUE
         )
##    print(result[[k]]$cpu.used)
    print(result.[[k]]$cpu.used)
}
print(rbind(result.[[1]]$cpu.used,
            result.[[2]]$cpu.used))
