## Export: inla.lattice2node.mapping inla.node2lattice.mapping inla.lattice2node inla.node2lattice inla.matrix2vector inla.vector2matrix

##!\name{lattice2node}
##!\alias{lattice2node}
##!\alias{inla.lattice2node}
##!\alias{node2lattice}
##!\alias{inla.node2lattice}
##!\alias{node2lattice.mapping}
##!\alias{inla.node2lattice.mapping}
##!\alias{lattice2node.mapping}
##!\alias{inla.lattice2node.mapping}
##!\alias{matrix2vector}
##!\alias{vector2matrix}
##!\alias{inla.matrix2vector}
##!\alias{inla.vector2matrix}
##!
##!\title{Functions to define mapping between a lattice and nodes}
##!
##!\description{These functions define mapping in between two-dimensional
##!  indices on a lattice and the one-dimensional node representation used
##!  in \code{inla}.
##!
##!  The mapping from node to lattice follows the default \code{R}
##!  behaviour (which is column based storage), and \code{as.vector(A)} and
##!  \code{matrix(a, nrow, ncol)} can be used instead of
##!  \code{inla.matrix2vector} and \code{inla.vector2matrix}.
##!}
##!
##!\usage{
##!inla.lattice2node.mapping(nrow, ncol)
##!inla.node2lattice.mapping(nrow, ncol)
##!inla.lattice2node(irow, icol, nrow, ncol)
##!inla.node2lattice(node, nrow, ncol)
##!inla.matrix2vector(a.matrix)
##!inla.vector2matrix(a.vector, nrow, ncol)
##!}
##!
##!\arguments{
##!
##!  \item{nrow}{Number of rows in the lattice.}
##!
##!  \item{ncol}{Number of columns in the lattice.}
##!
##!  \item{irow}{Lattice row index, between \code{1} and \code{nrow}}
##!
##!  \item{icol}{Lattice column index, between \code{1} and \code{ncol}}
##!
##!  \item{node}{The node index, between \code{1} and \code{ncol*nrow}}
##!
##!  \item{a.matrix}{is a matrix to be mapped to a vector using internal
##!  representation defined by \code{inla.lattice2node}}
##!
##!  \item{a.vector}{is a vector to be mapped into a matrix using the
##!    internal representation defined by \code{inla.node2lattice}}
##!}
##!
##!\value{\code{inla.lattice2node.mapping} returns the hole mapping as a
##!  matrix, and \code{inla.node2lattice.mapping} returns the hole mapping
##!  as \code{list(irow=..., icol=...)}. \code{inla.lattice2node} and
##!  \code{inla.node2lattice} provide the mapping for a given set of
##!  lattice indices and nodes. \code{inla.matrix2vector} provide the
##!  mapped vector from a matrix, and \code{inla.vector2matrix} provide the
##!  inverse mapped matrix from vector.}
##!  
##!%%
##!
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!
##!\seealso{\link{inla}}
##!
##!\examples{
##!## write out the mapping using the two alternatives
##!nrow = 2
##!ncol = 3
##!mapping = inla.lattice2node.mapping(nrow,ncol)
##!
##!for (i in 1:nrow){
##!    for(j in 1:ncol){
##!        print(paste("Alt.1: lattice index [", i,",", j,"] corresponds",
##!                    "to node [", mapping[i,j],"]", sep=""))
##!    }
##!}
##!
##!for (i in 1:nrow){
##!    for(j in 1:ncol){
##!        print(paste("Alt.2: lattice index [", i,",", j,"] corresponds to node [",
##!                    inla.lattice2node(i,j,nrow,ncol), "]", sep=""))
##!    }
##!}
##!
##!inv.mapping = inla.node2lattice.mapping(nrow,ncol)
##!for(node in 1:(nrow*ncol))
##!   print(paste("Alt.1: node [", node, "] corresponds to lattice index [",
##!               inv.mapping$irow[node], ",",
##!               inv.mapping$icol[node],"]", sep=""))
##!
##!for(node in 1:(nrow*ncol))
##!   print(paste("Alt.2: node [", node, "] corresponds to lattice index [",
##!               inla.node2lattice(node,nrow,ncol)$irow[1], ",",
##!               inla.node2lattice(node,nrow,ncol)$icol[1],"]", sep=""))
##!
##!## apply the mapping from matrix to vector and back
##!n = nrow*ncol
##!z = matrix(1:n,nrow,ncol)
##!z.vector = inla.matrix2vector(z)  # as.vector(z) could also be used
##!print(mapping)
##!print(z)
##!print(z.vector)
##!
##!## the vector2matrix is the inverse, and should give us the z-matrix
##!## back. matrix(z.vector, nrow, ncol) could also be used here.
##!z.matrix = inla.vector2matrix(z.vector, nrow, ncol) 
##!print(z.matrix)
##!}

`inla.lattice2node.mapping` = function(nrow, ncol)
{
    ## return a matrix with the mapping

    stopifnot( nrow > 0 && ncol > 0 )

    mapping = matrix(NA, nrow=nrow, ncol=ncol)
    for(i in 1:nrow) {
        j = 1:ncol
        mapping[i, j] = inla.lattice2node(i, j, nrow, ncol)
    }
    return (mapping)
}

`inla.node2lattice.mapping` = function(nrow, ncol)
{
    stopifnot( nrow > 0 && ncol > 0 )

    return (inla.node2lattice(1:(nrow*ncol), nrow, ncol))
}

`inla.lattice2node` = function(irow, icol, nrow, ncol)
{
    ## convert from a lattice point (irow, icol) to a node-number in
    ## the graph; similar to the GMRFLib_lattice2node()-function in
    ## GMRFLib.  Indices goes from irow=1...nrow, to icol=1...ncol and
    ## node=1.....nrow*ncol.

    stopifnot( nrow > 0 && ncol > 0 )

    if (length(irow) == length(icol) && length(irow) > 1) {
        ## this makes it kind of 'vectorize' for two arguments...
        n = length(irow)
        k = 1:n
        return (sapply(k,
                       function(irow, icol, nrow, ncol, k) {
                           return (inla.lattice2node(irow[k], icol[k], nrow, ncol))
                       }, irow = irow, icol = icol, nrow = nrow, ncol = ncol))
    } else {
        return ((irow-1) + (icol-1)*nrow + 1)
    }
}

`inla.node2lattice` = function(node, nrow, ncol)
{
    ## convert from a node-number in the graph to a lattice point
    ## (irow, icol); similar to the GMRFLib_node2lattice()-function in
    ## GMRFLib.  Indices goes from irow=1...nrow, to icol=1...ncol and
    ## node=1.....nrow*ncol.

    stopifnot( nrow > 0 && ncol > 0 )

    icol = (node - 1) %/% nrow
    irow = (node - 1) %% nrow

    return (list(irow = irow+1, icol = icol+1))
}

`inla.matrix2vector` = function(a.matrix)
{
    ## utility function for mapping a matrix to inla's internal `node'
    ## representation by inla.lattice2node() and inla.node2lattice()

    if (!is.matrix(a.matrix))
        stop("Argument must be a matrix")

    return (as.vector(a.matrix))
}
`inla.vector2matrix` = function(a.vector, nrow, ncol)
{
    ## utility function for mapping from inla's internal `node'
    ## representation, inla.lattice2node() and inla.node2lattice(),
    ## and to a matrix

    n = length(a.vector)
    if (missing(nrow) && !missing(ncol))
        nrow = n %/% ncol
    if (!missing(nrow) && missing(ncol))
        ncol = n %/% nrow
    if (n != nrow*ncol)
        stop(paste("Length of vector", n, "does not equal to nrow*ncol", nrow*ncol))

    return (matrix(a.vector, nrow, ncol))
}
