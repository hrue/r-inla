inla.simplemvspde.create = function(mesh) 
{
	inla.require.inherits(mesh,"inla.mesh","mesh")
	
	simplemvspde.prefix = inla.fmesher.make.prefix(NULL,NULL)
	
	simplemvspde = (list(mesh=mesh,
		        n.spde = mesh$n, #number of dofs in each spde
			n.components = 2, #always!!
			n.theta = 6,  ##GENERALISE
			param.inla=list(),
		        f= list()
			))

	class(simplemvspde) = c("inla.simplemvspde","inla.model.class")

	fem = inla.fmesher.smorg(mesh$loc,
					mesh$graph$tv,
                	               	fem=2,
                         		output=list("c0", "g1", "g2"))
	M0 = fem$c0
	M1 = fem$g1
	M2 = fem$g2  #these can be generalised slightly to non-white rhs noise

	fmesher.write(inla.affirm.double(M0), simplemvspde.prefix, "M0")
    	fmesher.write(inla.affirm.double(M1), simplemvspde.prefix, "M1")
    	fmesher.write(inla.affirm.double(M2), simplemvspde.prefix, "M2")


	#set some  initial thetas
	mesh.range = (max(c(diff(range(mesh$loc[,1])),
                            diff(range(mesh$loc[,2])),
                            diff(range(mesh$loc[,3]))
                            )))
	kappa0 = sqrt(8)/(mesh.range*0.2)
        tau0 = 1/sqrt(4*pi*kappa0^2)/1.0
	
	theta.initial=c(-2*log(kappa0),2*log(kappa0),2*log(kappa0),2*log(tau0),0.0,2*log(tau0))  ## set b21 = 0 initially!


	param.inla = list( n = mesh$n, n.theta=6,
		     	   M0=M0,M1=M1,M2=M2,
			   transform="identity",
			   theta.initial=theta.initial,
			   fixed=FALSE,
			   theta.fixed=NULL)
	simplemvspde$param.inla=param.inla

	simplemvspde$f = ( list( model="simplemvspde",
		       	   	 simplemvspde.prefix=simplemvspde.prefix,
				 simplemvspde.transform = "identity",
				 hyper.default=list(
					theta1=list(prior="normal",param=c(theta.initial[1],0.1)),
					theta2=list(prior="normal",param=c(theta.initial[2],0.1)),
					theta3=list(prior="normal",param=c(theta.initial[3],0.1)),
					theta4=list(prior="normal",param=c(theta.initial[4],0.1)),
					theta5=list(prior="normal",param=c(theta.initial[5],0.1)),
					theta6=list(prior="normal",param=c(theta.initial[6],0.1))
				 )
	))


	return(invisible(simplemvspde))
}

inla.simplemvspde.precision = function(simplemvspde,kappas,bs)
{
	 
	#precision =NULL
	M0 = simplemvspde$param.inla$M0 #c0
	M1 = simplemvspde$param.inla$M1 #g1
	
	n = dim(M0)[1]
	print(n)
	Zero = sparseMatrix(i=1,j=1,x=0,dims=c(n,n))
	L = rBind(
	  cBind(bs[1]*(kappas[1]*M0 + M1), Zero),
	  cBind(bs[2]*(kappas[2]*M0 + M1),bs[3]*( kappas[3]*M0 + M1))
	)
	Q = sparseMatrix(i=1:(2*n),j=1:(2*n),x=c(1/diag(M0),1/diag(M0)),dims=c(2*n,2*n))
	precision = t(L)%*%Q%*%L
	precision = 0.5*(precision + t(precision))
 	return(invisible(precision))
}