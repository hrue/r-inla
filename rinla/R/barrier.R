## Export: inla.barrier.pcmatern inla.barrier.polygon

##! \name{inla.barrier}
##! \alias{inla.barrier}
##! \alias{barrier}
##! \alias{inla.barrier.pcmatern}
##! \alias{barrier.pcmatern}
##! \alias{inla.barrier.polygon}
##! \alias{barrier.polygon}
##! \alias{inla.barrier.q}
##! \alias{barrier.q}
##! \alias{inla.barrier.fem}
##! \alias{barrier.fem}
##! 
##! \title{Functions for defining the Barrier models}
##! 
##! \description{Functions for defining Barrier models as an \code{inla rgeneric} model}
##!
##! \usage{
##! inla.barrier.pcmatern(mesh, barrier.triangles, prior.range,
##!                       prior.sigma, range.fraction=0.2)
##! inla.barrier.polygon(mesh, barrier.triangles, Omega=NULL)
##! inla.barrier.q(fem, ranges, sigma=1)
##! inla.barrier.fem(mesh, barrier.triangles, Omega=NULL)
##! }
##! \arguments{
##!   \item{mesh}{The mesh to build the model on, from inla.mesh.2d}
##!   \item{barrier.triangles}{The numerical ids of the triangles that make up the barrier area}
##!   \item{prior.range}{2 parameters \code{(range0,Prange)} for the prior spatial range. 
##!         If \code{Prange} is \code{NA}, then \code{range0} is used as a fixed range value (TODO).}
##!   \item{prior.sigma}{2 parameters \code{(sig0,Psig)} for the prior marginal standard deviation sigma. 
##!         If \code{Psig} is \code{NA}, then \code{sig0} is used as a fixed sigma value (TODO).}
##!   \item{range.fraction}{The length of the spatial range inside the barrier area,
##!                         as a fraction of the range parameter.}
##!   \item{Omega}{Advanced option for creating a set of permeable barriers (not documented)}
##! }
##! \details{
##!     This model is described in the ArXiv preprint arXiv:1608.03787.
##!     For examples, see \url{https://haakonbakka.bitbucket.io/btopic107.html}.
##! }
##!\value{%%
##!  \code{inla.barrier.pcmatern} gives the (rgeneric) model object for fitting the model in INLA,
##!  \code{inla.barrier.polygon} gives the polygon around the barrier (mainly for plotting),
##!  \code{inla.barrier.q} is an internal method producing the Q matrix from a result of inla.barrier.fem,
##!  \code{inla.barrier.fem} is an internal method producing the Finite Element matrices.
##! }
##! \seealso{inla.spde2.pcmatern}
##! \author{Haakon Bakka \email{bakka@r-inla.org}}


inla.barrier.pcmatern = function(mesh, barrier.triangles, prior.range, prior.sigma, range.fraction=0.2)
{
    ## Give default values if absolutely needed
    if (missing(prior.range)) {
        warning("Arbitrary prior values chosen automatically. This may suffice for a first attempt, 
            but should be changed in any serious analysis.")
        prior.range = c(diff(range(mesh$loc[ ,1]))/5, 0.5)
    }
    
    if (missing(prior.sigma)) {
        warning("Arbitrary prior values chosen automatically. This may suffice for a first attempt, 
            but should be changed in any serious analysis.")
        prior.sigma = c(1, 0.5)
    }
    
    ## INPUT verification ###
    stopifnot(class(mesh) == 'inla.mesh')
    stopifnot(range.fraction > 0.000001)

    ## FUNCTIONS FOR RGENERIC MODEL SETUP ###
    
    dt.create.prior.log.exp = function (prior.param) {
        ## This is the log-prior for the internal parametrisation theta
        ## Both log of probability and log of exponential dist
        ## theta = log(sigma), log(range1), log(range2), ...
        ## Input
        ## parameters are the lambdas in the exponential distribution
        ## - the first is for sigma
        ## - the second is for all the ranges
        
        ## Move to current scope (environment)
        prior.param = prior.param
        
        log.prior = function(theta) {
            lambda0 = prior.param[1]
            lambda1 = prior.param[2]
            ntheta = length(theta)
            
            ## Prior for standard deviation
            val = 0 + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
            
            ## Prior for range(s)
            for (i in 2:ntheta) {
                val = val + log(lambda1) - lambda1*exp(-theta[i]) + -theta[i]
            }
            return(val)
        }
        ## - this environment includes the prior parameters
        return(log.prior)
    }
    
    
    ## - this function is the model component definition in the rgeneric inla framework
    ## - see inla.doc('rgeneric')
    barrier.rgeneric.model = function(cmd = c("graph", "Q", "mu", "initial",
                                              "log.norm.const", "log.prior", "quit"),
                                      theta = NULL) {
        ## Input
        ## theta
        
        ## Assumed functions
        ## log.prior(theta)
        
        graph = function(theta) {
            require(methods)
            ntheta = length(initial())
            G1 = Q(theta=(1:ntheta)/3.217233456)
            G1[G1 != 0] = 1
            G2 = Q(theta=(1:ntheta)^2/12.1543534)
            G2[G2 != 0] = 1
            
            return (G1+G2)
        }

        Q = function(theta) {
            return(inla.barrier.q(fem=fem, ranges=exp(theta[2])*c(1, range.fraction), sigma=exp(theta[1])))
        }
        
        mu = function(theta) numeric(0)
        
        initial = function(theta) {
            initial.theta = rep(0, 2)
            return (initial.theta)
        }
        
        log.norm.const = function(theta) numeric(0)
        
        quit <- function(theta) invisible()
        
        val = do.call(match.arg(cmd), args = list(theta))
        return (val) 
    }
    
    ## Create a valid Omega from the barrier.triangles
    barrier.triangles = unique(barrier.triangles)
    t = length(mesh$graph$tv[,1])

    ## - all triangles
    remaining = setdiff(1:t, barrier.triangles)
    Omega = list(remaining, barrier.triangles)
    
    ## Create the barrier model
    if (!is.na(prior.sigma[2]) && !is.na(prior.range[2])) {
        log.prior = dt.create.prior.log.exp(
            prior.param = c(-log(prior.sigma[2])/prior.sigma[1], -log(prior.range[2])*prior.range[1]))
        ## - The prior parameters are the lambdas in the exponential 
        ##   priors for standard deviation and inverse-range
        ## - the first is log(prob)/exceed, the second log(prob)*exceed
        ## - the second is exponential for inverse range, therefore multiplication!
        fem = inla.barrier.fem(mesh, barrier.triangles = barrier.triangles)
        barrier.model = inla.rgeneric.define(model = barrier.rgeneric.model,
                                             log.prior=log.prior, inla.barrier.q=inla.barrier.q, 
                                             fem=fem, range.fraction = range.fraction)
    } else if (!is.na(prior.sigma[2]) && is.na(prior.range[2])) {
        stop("Input not supported (TODO)")
    } else {
        stop("Input not supported (TODO)")
    }

    return(barrier.model)
}

inla.barrier.polygon = function(mesh, barrier.triangles, Omega=NULL)
{
    ## - constructs SpatialPolygons for the different subdomains (areas)
    stopifnot(class(mesh) == 'inla.mesh')
    ## - requires an inla mesh to work
    library(rgeos)
    ## - required package
    
    if (missing(barrier.triangles)) { 
        ## Use Omega
    } else {
        ## Create a valid Omega
        barrier.triangles = unique(barrier.triangles)
        t = length(mesh$graph$tv[,1])
        ## - all triangles
        remaining = setdiff(1:t, barrier.triangles)
        if (!is.null(Omega)) warning("Omega is replaced by barrier.triangles")
        Omega = list(remaining, barrier.triangles)
    }
    
    Omega.SP.list = list()
    for (j in 1:length(Omega)) {
        poly.list = list()
        for (tri in Omega[[j]]) {
            px = mesh$graph$tv[tri, ]
            temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            poly.list = c(poly.list , Polygon(rbind(temp[ ,1:2], temp[1, 1:2]), hole=F))
        }
        mesh.polys = SpatialPolygons(list(Polygons(poly.list, ID='noid')))
        Omega.SP.list[[j]] = gUnaryUnion(mesh.polys)
    }
    
    if (missing(barrier.triangles)) {
        return(Omega.SP.list)
    } else {
        return(Omega.SP.list[[2]])
    }
}

### PRECISION MATRIX FUNCTIONS ###
## I.e. Solve the differential equation (the SPDE)

inla.barrier.q <- function(fem, ranges, sigma=1)
{
    ## - This function computes a specific precision matrix
    ## - fem represents the Barrier model or the Different Terrains model
    ## - fem contains all the needed matrices to solve the SPDE
    ## - the ranges and sigma are the hyperparameters that determine Q
    
    if (is.null(ranges)) stop("ranges cannot be NULL")
    if (any(is.na(ranges))) stop("No range can be NA")
    
    xi = length(ranges)
    if (xi != length(fem$D)) {
        print('dt.precision has encountered an error. Will stop.')
        stop ('Ranges do no correspond to fem')
    }
    if (any(ranges < 0.001)) {
        warning('This hyper parameter value will probably fail. A very small maximum edge length needed in the mesh.')
    }
    
    Cdiag = ranges[1]^2* fem$C[[1]] # already raised to power 2
    if (xi > 1) {
        for (k in 2:xi) {
            Cdiag = Cdiag + ranges[k]^2*fem$C[[k]]
        }
    }
    N = length(Cdiag)
    Cinv = sparseMatrix(i=1:N, j = 1:N, x=1/Cdiag, dims = c(N,N), giveCsparse=FALSE) 
    
    A = fem$I  
    for (k in 1:xi) {
        A = A + (ranges[k]^2/8)*fem$D[[k]]
    }
    
    Q = t(A)%*%Cinv%*%A*(1/sigma^2)/(pi/2) * 4 
    Q = inla.as.dgTMatrix(Q)
    return (Q)
}


`inla.barrier.fem` = function(mesh, barrier.triangles, Omega = NULL)
{
    stopifnot(class(mesh) == 'inla.mesh')
    
    dt.fem.matrices <- function(mesh, Omega) {
        ## - This function computes the Finite Element matrices
        ## - - this is needed to compute the precision matrix Q later
        
        xi = length(Omega)
        fem = list()
        fem$I = dt.fem.identity(mesh)
        fem$D = list()
        fem$C = list()
        for (k in 1:xi) {
            fem$D[[k]] = dt.fem.laplace(mesh, Omega[[k]])
        }
        for (k in 1:xi) {
            fem$C[[k]] = dt.fem.white(mesh, Omega[[k]])
        }
        fem$hdim = xi
        return(fem)
    }
    
    dt.fem.white <- function(mesh, subdomain) {
        ## - This function computes the Finite Element matrix of the white noise on ONE subdomain (area)
        ## - This matrix is a diagonal matrix
        
        ## Pre-allocation
        Ck = rep(0, mesh$n)
        
        for (t in subdomain) {
            ## Node indexes for triangle t:
            px = mesh$graph$tv[t, ]
            temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            p1 = t(t(temp[1, c(1,2)]))
            p2 = t(t(temp[2, c(1,2)]))
            p3 = t(t(temp[3, c(1,2)]))
            Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
            area=  abs(det(Ts)) * 0.5
            for (i in 1:3) {
                Ck[px[i]] = Ck[px[i]] + area
            }
        }
        return (Ck)
    }
    
    dt.fem.identity <- function(mesh) {
        ## - this function computes the Finite Element matrix for the '1' in the SPDE (the identity operator)
        ## - this operator does not depend on the subdomains
        
        ## Preallocation
        len = length(mesh$graph$tv[,1])
        index.i = rep(0,len * 6)
        index.j = rep(0,len * 6)
        Aij = rep(0,len * 6) # matrix values
        counter = 1
        
        for (t in 1:len) {
            ## Node indexes for triangle t:
            px = mesh$graph$tv[t, ]
            temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            p1 = t(t(temp[1, c(1,2)]))
            p2 = t(t(temp[2, c(1,2)]))
            p3 = t(t(temp[3, c(1,2)]))
            Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
            twiceArea=  abs(det(Ts))
            
            for (i in 1:3) {
                index.i[counter] = px[i]
                index.j[counter] = px[i] # same
                Aij[counter] = (twiceArea)*1/12
                counter = counter + 1
            }
            
            for (i in 1:2) {
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

    dt.fem.laplace <- function(mesh, subdomain) {
        ## - This function computes the Finite Element matrix of the laplace operator on ONE subdomain (area)
        ## - This matrix is very sparse
        ## The nabla phi's
        Nphix = rbind(c(-1,-1), c(1,0), c(0,1))
        len = length(subdomain)
        index.i = rep(0,len * 9)
        index.j = rep(0,len * 9)
        Aij = rep(0,len * 9) # matrix values
        counter = 1
        
        for (tri in subdomain) {
            px = mesh$graph$tv[tri, ]
            temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            p1 = t(t(temp[1, c(1,2)]))
            p2 = t(t(temp[2, c(1,2)]))
            p3 = t(t(temp[3, c(1,2)]))
            Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
            TTTinv = solve(t(Ts)%*%Ts)
            area=  abs(det(Ts)) * 0.5
            
            for (k in 1:3) {
                for (m in 1:3) {
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
    
    if (missing(barrier.triangles) && is.null(Omega)) stop("Input barrier triangles")
    
    if (missing(barrier.triangles)) { 
        ## Use Omega
    } else {
        ## Create a valid Omega
        barrier.triangles = unique(barrier.triangles)
        t = length(mesh$graph$tv[,1])
        ## - all triangles
        remaining = setdiff(1:t, barrier.triangles)
        if (!is.null(Omega)) warning("Omega is replaced by barrier.triangles")
        Omega = list(remaining, barrier.triangles)
    }
    
    fem = dt.fem.matrices(mesh, Omega)
    return(fem)
}

