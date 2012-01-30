#require(INLA)
require(rgl) #For plotting!

require(VGAM) #for rfrechet

# Make some fake data
N=500

Ntrials = rep(100,N)
points = matrix(runif(2*N),ncol=2)

#The 'true' surface

fcn = function(loc){
	return( cos(loc[,1])*sin(4*loc[,2]))
	}
	
noisy_data = rfrechet(N, location=fcn(points),shape=1.0)

## Build a mesh

bnd = inla.mesh.segment(matrix(c(0,0,1,0,1,1,0,1),ncol=2,byrow=TRUE))

mesh = inla.mesh.create(points,boundary=bnd,refine=list(max.edge=0.1))

plot(mesh)

##Make the SPDE

spde = inla.spde.create(mesh,model="imatern")

## Fit the model

fake_data = list(y=noisy_data, data_points = mesh$idx$loc)

hyper = list( prec = list(initial=1) ,  gev = list(intial=1))

formula = y ~ f(data_points, model=spde)-1

r = inla(formula, family="gev", data=fake_data,verbose=TRUE,control.mode = list(theta= c(1,1,1), restart=TRUE))

if(FALSE){
## And some output
summary(r)
plot(r, plot.random.effects=FALSE)

plot(old.mesh.class(mesh),r$summary.random$data_points$mean)

#Calculate the MSE

print(sqrt(var(r$summary.random$data_points$mean - fcn(mesh$loc))))
}
