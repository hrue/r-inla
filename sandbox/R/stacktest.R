n = 1000
max.edge = 0.05
mesh =
    inla.mesh.create(boundary=(inla.mesh.segment(matrix(c(0,1,1,0,0,0,1,1),
                                                        4,2), is.bnd=TRUE)),
                     refine=list(max.edge=max.edge))
spde =
    inla.spde2.matern(mesh)

loc = matrix(runif(n*2),n,2)
covar1 = matrix(runif(n),n,1)
covar2 = matrix(runif(n),n,1)
A = inla.spde.make.A(mesh, loc=loc)
y1 = covar1+rnorm(n)*0.01
y2 = covar2+rnorm(n)*0.01

spatial = inla.spde.make.index("spatial", n.mesh=mesh$n)

stack =
    inla.stack(data = list(y=y1),
               A = list(A, 1),
               effect =
               list(spatial,
                    list(covar=covar1)))


stack = list()
stack[[1]] =
    inla.stack(data = list(y=y1),
               A = list(A, 1),
               effect =
               list(spatial,
                    list(covar=covar1)))
stack[[2]] = inla.stack(stack[[1]], stack[[1]])
stack[[3]] =
    inla.stack(data = list(y=c(y1,y2)),
               A = list(rBind(A,A), Diagonal(n*2,1)),
               effect =
               list(inla.spde.make.index("spatial", n.field=mesh$n),
                    list(covar=c(covar1,covar2))))

stack[[4]] =
    inla.stack(data = list(y=y1),
               A = list(A, covar1),
               effect =
               list(inla.spde.make.index("spatial", n.field=mesh$n),
                    list(covar=1)))
stack[[5]] =
    inla.stack(data = list(y=c(y1,y2)),
               A = list(rBind(A,A), rBind(covar1,covar2)),
               effect =
               list(inla.spde.make.index("spatial", n.field=mesh$n),
                    list(covar=1)))

formula = y ~ -1 + covar + f(spatial, model=spde)

result = list()
for (k in 1:length(stack)) {
    result[[k]] =
        inla(
             formula, data=stack[[k]]$data, family="gaussian",
             control.predictor=list(A=stack[[k]]$A)
         )
    print(result[[k]]$cpu.used)
}
print(rbind(result[[1]]$cpu.used,
            result[[2]]$cpu.used,
            result[[3]]$cpu.used,
            result[[4]]$cpu.used,
            result[[5]]$cpu.used))
