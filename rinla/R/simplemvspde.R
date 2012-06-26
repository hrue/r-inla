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
	// make this function
	precision =NULL
	return(invisible(precision))
}