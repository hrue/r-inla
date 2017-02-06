### DESCRIPTION ###
# This is a minimal set of fuctions, for doing the DT model

### DESCRIPTION ###
# These functions accompanies the script 'MinimalExample'
# These functions are needed to run the Barrier model
# These functions can also be used to run the more general Different Terrains model

### LOAD LIBRARIES ###
library(INLA)
library(rgeos)
library(fields)
library(rgdal)

### GENERAL FUNCTIONS ###
print.v = function(x, ...) {
  # - Only print if the verbose variable is true (or not exists)
  if(!exists('verbose') || !is.logical(verbose) || verbose) {
    print(x, ...)
  }
}


bakka.square.polygon = function(xlim, ylim, ret.SP = F, id=runif(1)){
  # - ret.SP=T : Return a SpatialPolygons object
  # - ret.SP=F : Return a Polygon object
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, c(corner2[1], corner1[2]), corner1), hole = FALSE)
  if (ret.SP) {
    return(SpatialPolygons(list(Polygons(list(poly), ID = id))))
  } else {
    return(poly)
  }
}

inla.mesh.add.posTri <- function(mesh){
  # - Add on two attributes to the mesh object
  mesh$t = length(mesh$graph$tv[,1])
  mesh$posTri = matrix(0, mesh$t, 2)
  for (t in 1:mesh$t){
    temp = mesh$loc[mesh$graph$tv[t, ], ]
    mesh$posTri[t,] = colMeans(temp)[c(1,2)] 
  }
  return(mesh)
}

dt.Omega <- function(list_of_subdomains, mesh){
  # - list_of_subdomains must be a list of numeric 'c(...)'
  # - This fun creates a legal Omega from just a list of subdomains
  # - It does not make sure every triangle is included somewhere
  # - But it does make sure there are no multiplicities
  
  if (length(list_of_subdomains[[1]])==0){
    stop('No proper first subdomain')
  }
  Omega = list()
  Omega[[1]] = list_of_subdomains[[1]]
  if (length(list_of_subdomains)==1){
    remainder = setdiff(1:mesh$t, Omega[[1]])
    if (length(remainder)>0){
      Omega[[length(Omega)+1]] = remainder
    }
    return (Omega)
  }
  usedTriangles = Omega[[1]]
  for (i in 2:length(list_of_subdomains)){
    Omega[[i]] = setdiff(list_of_subdomains[[i]], usedTriangles)
    usedTriangles = union(usedTriangles, Omega[[i]])
  }
  ## Remove empty Omega parts
  for (i in 1:length(list_of_subdomains)){
    if (length(Omega[[length(Omega)]])==0) Omega[[length(Omega)]] = NULL
  }
  
  ## Add on remainder
  remainder = setdiff(1:mesh$t, usedTriangles)
  if (length(remainder)>0){
    Omega[[length(Omega)+1]] = remainder
  }
  return (Omega)
}


dt.polygon.omega = function (mesh, Omega) {
  # - constructs SpatialPolygons for the different subdomains (areas)
  stopifnot(class(mesh) == 'inla.mesh')
  # - requires an inla mesh to work
  
  Omega.SP.list = list()
  for (j in 1:length(Omega)) {
    poly.list = list()
    for (tri in Omega[[j]]){
      px = mesh$graph$tv[tri, ]
      temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
      poly.list = c(poly.list , Polygon(rbind(temp[ ,1:2], temp[1, 1:2]), hole=F))
    }
    mesh.polys = SpatialPolygons(list(Polygons(poly.list, ID='noid')))
    Omega.SP.list[[j]] = gUnaryUnion(mesh.polys)
  }
  return(Omega.SP.list)
}


dt.nearest.mesh.index = function (location, mesh, verbose = T) {
  # Uses the A-construction and picks one of the nearest mesh nodes
  
  A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(c(location[1], location[2]),1,2))
  node.id = which.max(A.tmp[1, ])
  
  return(node.id)
}



### PRECISION MATRIX FUNCTIONS ###
# I.e. Solve the differential equation (the SPDE)

dt.precision <- function(spde, ranges, sigma=1){
  # - This function computes a specific precision matrix
  # - spde is the Different Terrains model
  # - spde contains all the needed matrices to solve the SPDE
  # - the ranges and sigma are the hyperparameters that determine Q
  
  xi = length(ranges)
  if (xi != length(spde$D)){
    print('dt.precision has encountered an error. Will stop.')
    stop ('Ranges do no correspond to spde')
  }
  if (any(ranges < 0.001)){
    warning('This hyper parameter value will probably fail. A very small maximum edge length needed in the mesh.')
  }
  
  Cdiag = ranges[1]^2* spde$C[[1]] # already raised to power 2
  if (xi > 1){
    for (k in 2:xi){
      Cdiag = Cdiag + ranges[k]^2*spde$C[[k]]
    }
  }
  N = length(Cdiag)
  Cinv = sparseMatrix(i=1:N, j = 1:N, x=1/Cdiag, dims = c(N,N), giveCsparse=FALSE) 
  
  A = spde$I  
  for (k in 1:xi){
    A = A+ (ranges[k]^2/8)*spde$D[[k]]
  }
  
  Q = t(A)%*%Cinv%*%A*(1/sigma^2)/(pi/2) *4 
  Q = inla.as.dgTMatrix(Q)
  return (Q)
}


dt.FEMmatrices <- function(mesh, Omega){
  # - This function computes the Finite Element matrices
  # - - this is needed to compute the precision matrix Q later
  
  xi = length(Omega)
  spde = list()
  spde$I = dt.FEMidentity(mesh)
  spde$D = list()
  spde$C = list()
  for (k in 1:xi) {
    spde$D[[k]] = dt.FEMlaplace(mesh, Omega[[k]])
  }
  for (k in 1:xi){
    spde$C[[k]] = dt.FEMwhite(mesh, Omega[[k]])
  }
  spde$hdim = xi
  return(spde)
}


dt.FEMwhite <- function(mesh, subdomain){
  # - This function computes the Finite Element matrix of the white noise on ONE subdomain (area)
  # - This matrix is a diagonal matrix
  
  ## Pre-allocation
  Ck = rep(0, mesh$n)
  
  for (t in subdomain){
    ## Node indexes for triangle t:
    px = mesh$graph$tv[t, ]
    temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
    p1 = t(t(temp[1, c(1,2)]))
    p2 = t(t(temp[2, c(1,2)]))
    p3 = t(t(temp[3, c(1,2)]))
    Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
    area=  abs(det(Ts)) * 0.5
    for (i in 1:3){
      Ck[px[i]] = Ck[px[i]] + area
    }
  }
  return (Ck)
}


dt.FEMidentity <- function(mesh){
  # - this function computes the Finite Element matrix for the '1' in the SPDE (the identity operator)
  # - this operator does not depend on the subdomains
  
  ## Preallocation
  len = length(mesh$graph$tv[,1])
  index.i = rep(0,len * 6)
  index.j = rep(0,len * 6)
  Aij = rep(0,len * 6) # matrix values
  counter = 1
  
  for (t in 1:len){
    ## Node indexes for triangle t:
    px = mesh$graph$tv[t, ]
    temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
    p1 = t(t(temp[1, c(1,2)]))
    p2 = t(t(temp[2, c(1,2)]))
    p3 = t(t(temp[3, c(1,2)]))
    Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
    twiceArea=  abs(det(Ts))
    
    for (i in 1:3){
      index.i[counter] = px[i]
      index.j[counter] = px[i] # same
      Aij[counter] = (twiceArea)*1/12
      counter = counter + 1
    }
    
    for (i in 1:2){
      for (j in (i+1):3)
        index.i[counter] = px[i]
      index.j[counter] = px[j]
      Aij[counter]= (twiceArea)*1/24
      counter=counter + 1
      # symmetry:
      index.i[counter] = px[j]
      index.j[counter] = px[i]
      Aij[counter]= (twiceArea)*1/24
      counter=counter + 1
    }
  }
  
  I = sparseMatrix(i=index.i, j = index.j, x=Aij, dims = c(mesh$n, mesh$n), giveCsparse=FALSE) 
  
  return (I)
}



dt.FEMlaplace <- function(mesh, subdomain){
  # - This function computes the Finite Element matrix of the laplace operator on ONE subdomain (area)
  # - This matrix is very sparse
  
  # The nabla phi's
  Nphix = rbind(c(-1,-1), c(1,0), c(0,1))
  len = length(subdomain)
  index.i = rep(0,len * 9)
  index.j = rep(0,len * 9)
  Aij = rep(0,len * 9) # matrix values
  counter = 1
  
  for (tri in subdomain){
    px = mesh$graph$tv[tri, ]
    temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
    p1 = t(t(temp[1, c(1,2)]))
    p2 = t(t(temp[2, c(1,2)]))
    p3 = t(t(temp[3, c(1,2)]))
    Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
    TTTinv = solve(t(Ts)%*%Ts)
    area=  abs(det(Ts)) * 0.5
    
    for (k in 1:3){
      for (m in 1:3){
        tmp = (3*m+k-4)*length(subdomain)
        index.i[(tmp + counter)] = px[k]
        index.j[(tmp + counter)] = px[m]
        Aij[(tmp + counter)] = area*Nphix[k, c(1,2) ]%*%TTTinv%*%as.matrix(Nphix[m, c(1,2) ])
      }
    }
    counter = counter + 1
  }
  
  Dk = sparseMatrix(i=index.i, j = index.j, x=Aij, dims = c(mesh$n, mesh$n), giveCsparse=FALSE) 
  return (Dk)
}



### FUNCTIONS FOR RGENERIC MODEL SETUP ###

## Check
if (!exists('dt.precision')) {
  stop('Error: You need additional functions')
}


dt.rgeneric.dt.model = function(
  # - this function is the model component definition in the rgeneric inla framework
  # - see demo(rgeneric)
  # - see inla.doc('rgeneric')
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL,
  args = NULL)
  
{
  
  interpret.theta = function(n, ntheta, theta) # not used right now
  {
    return (list(sigma = exp(theta[1L]),
                 r = exp(theta[2:ntheta])))
  }
  
  Q = function(n, ntheta, theta)
  {
    ## returns the precision matrix for given parameters
    
    t = interpret.theta(n, ntheta, theta) # parameters external scale
    
    Q = DT$Q(ranges = t$r, sigma=t$sigma) 
    
    Q = inla.as.dgTMatrix(Q)
    return (Q)
  }
  
  
  graph = function(n, ntheta , theta ) 
  {
    require(methods)
    
    G = Q(n, ntheta, (1:ntheta)/3.217233456)
    # - this should never give zeroes by accident
    # - the number 3.21... is an arbitrary choice
    
    G[G != 0] = 1
    # - not needed, as the values are interpreted as 1 or 0 anyway
    return (G)
  }
  
  
  log.norm.const = function(n, ntheta, theta)
  {
    return(numeric(0))
  }
  
  log.prior = function(n, ntheta, theta)
  {
      return(DT$prior(theta))
  }
  
  initial = function(n, ntheta, theta)
  {
    if(!is.null(DT$initial)) {
      stopifnot(length(DT$initial)==ntheta)
      return(DT$initial)
    }
    return (init = rep(1, ntheta)) 
  }
  
  quit = function(n, ntheta, theta)
  {
    return (invisible())
  }
  
  mu = function(n, ntheta, theta) 
  {
    return (rep(0, n))
  }
  
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list(
    n = as.integer(args$n),
    ntheta = as.integer(args$ntheta), 
    theta = theta))
  
  return (val)  
} # end function definition




















dt.rgeneric.environment = function (spde.dt, initial.values = NULL, fixed.ranges = NULL, copy.ranges = NULL, 
                                    verbose=F, prior.function = NULL, prior.parameters = c(1,1))
{
  # - This function creates the environment for the rgeneric definition of the Difficult Terrain model component
  
  ## Input
  # fixed.ranges
  # - NULL or c(T,F) list of values for if the hyper shall be fixed to a value
  # initial.values
  # - Initial values for setting the hypers 
  # - Internal scale: log(sigma) log(ranges)
  # - will be overridden by control.mode in the inla(...) statement
  # prior.function
  # - ust give the log prior as a function of the (internal) hyper-parameters
  
  # Not implemented: copy.ranges
  # - NULL or c(values) list of values smaller than index for an index to copy
  # - The purpose is to be able to use this to run a stationary model without creating another spde.dt variable
  if (!is.null(copy.ranges)) stop('Not implemented')
  
  DT = list()
  # - When inla(...) starts up the rgeneric R-process, this object is loaded
  # - This contains the necessary input to define the DT.rgeneric.dt.model
  # - This is the different terrains environment
  
  DT$spde = spde.dt
  DT$dt.precision = dt.precision
  if (is.null(fixed.ranges)){
    DT$Q  = function(ranges, sigma=1){
      return(DT$dt.precision(DT$spde, ranges, sigma))
    }
    
    ntheta = 1+length(spde.dt$D)
    # - create variable needed for rgeneric
    # - ntheta means number of thetas
  } else {
    # we have to fix some hyperparameters
    stopifnot(length(fixed.ranges) == length(spde.dt$D))
    DT$Q  = function(ranges, sigma=1){
      stopifnot(length(ranges)==sum(is.na(fixed.ranges)))
      internal.ranges = fixed.ranges
      internal.ranges[is.na(internal.ranges)] = ranges
      return(DT$dt.precision(DT$spde, internal.ranges, sigma))
    }
    ntheta = 1+sum(is.na(fixed.ranges))
    # - create variable needed for rgeneric
    # - ntheta means number of thetas
  }
  
  
  
  DT$prior.parameters = prior.parameters
  
  if (!is.null(prior.function)){ # the user specified function
    DT$prior = prior.function
  } else { # use the default prior function
    DT$prior = function(theta) {
      lambda0 = DT$prior.parameters[1]
      lambda1 = DT$prior.parameters[2]
      
      ## Prior for standard deviation
      val = 0 + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
      
      ## Prior for range(s)
      for (i in 2:ntheta) {
        val = val + log(lambda1) - lambda1*exp(-theta[i]) + -theta[i]
      }
      return(val)
    }
  }
  
  DT$verbose = verbose
  DT$ntheta = ntheta
  DT$initial = initial.values
  
  return (DT)
}









