require(INLA)
remove(list=ls())
inla.my.update()
mesh=
    inla.mesh.create.helper(matrix(c(0,1,1,0,0,0,1,0.99),4,2),
                            max.edge=c(0.02,0.5))
proj = inla.mesh.projector(mesh,dims=c(200,200))
for (al in seq(2,0.1,length.out=20)) {
    spde=inla.spde2.matern(mesh,alpha=al)
    Q=inla.spde2.precision(spde,theta=c(0,log(sqrt(8*max(0.5,al-1))/0.1)))
    field=inla.spde.sample(Q,seed=1L)
    image(inla.mesh.project(proj,field=field),
          main=paste("alpha =", al))
}


for (al in seq(2,0.1,length.out=20)) {
    nu=al-1
    spde=inla.spde2.matern(mesh,alpha=al)
    ka = sqrt(8*max(0.5,nu))/0.2
    s2 = inla.matern.cov(nu,ka,0,d=2)
    Q=inla.spde2.precision(spde,theta=c(0,log(ka)))
    A=inla.mesh.project(mesh,loc=matrix(c(0.5,0.5,0),1,3))$A
    S=rowSums(solve(Q, t(A)))

    H = sqrt((mesh$loc[,1]-0.5)^2+(mesh$loc[,2]-0.5)^2)
    h = (0:1000)/1000
    if (al>=1) {
        plot(h, inla.matern.cov(nu,ka,h,d=2), type="l",
             main=(paste("alpha =", al,
                         ", relative error =",
                         max(abs(S[H<0.4]-inla.matern.cov(nu,ka,H[H<0.4],d=2)))/s2)))
        points(H, S)
    } else {
        plot(h, inla.matern.cov(nu,ka,h,d=2), type="l",
             main=(paste("alpha =", al)))
        points(H, S)
    }

##    plot(H, abs(S/inla.matern.cov(nu,ka,H,d=2)-1), log="y")
}
