##!\name{graph2matrix}
##!\alias{graph2matrix}
##!\alias{inla.graph2matrix}
##!\alias{spy}
##!\alias{inla.spy}
##!\title{Construct a neighbour-matrix from a \code{graph}}
##!\description{Construct a neighbour-matrix from a \code{graph} and disaply it}
##!\usage{
##!Q = inla.graph2matrix(graph)
##!inla.spy(graph, reordering = NULL)
##!}
##!\arguments{
##!    \item{graph}{An \code{inla.graph}-object, a (sparse) symmetric matrix, a filename containing the graph,
##!                 or a list or collection of characters and/or numbers defining the graph.}
##!    \item{g}{An \code{inla.graph}-object}
##!    \item{Q}{An (possible) sparse symmtric matrix}
##!    \item{reordering}{A possible reordering. Typical the one obtained from a \code{inla}-call,  \code{result$misc$reordering}.}
##!}
##!\value{
##!     \code{inla.graph2matrix} returns a sparse symmetric matrix where the non-zero pattern is defined by the \code{graph}.
##!     The \code{inla.spy} function, plots a binary image of a \code{graph}. The \code{reordering} argument
##!     is typically the reordering used by \code{inla}, found in \code{result$misc$reordering}.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{
##!    \code{\link{inla.read.graph}}
##!}
##!\examples{
##!n = 50
##!Q = matrix(0, n, n)
##!idx = sample(1:n, 2*n, replace=TRUE)
##!Q[idx, idx] = 1
##!diag(Q) = 1
##!g = inla.read.graph(Q)
##!QQ = inla.graph2matrix(g)
##!inla.spy(QQ)
##!print(all.equal(as.matrix(Q), as.matrix(QQ)))
##!
##!g.file = inla.write.graph(g)
##!inla.dev.new()
##!inla.spy(g.file)
##!
##!g = inla.read.graph(g.file)
##!inla.dev.new()
##!inla.spy(g)
##!
##!inla.dev.new()
##!inla.spy(3, 1, "1 2 2 1 1 3 0")
##!inla.dev.new()
##!inla.spy(3, 1, "1 2 2 1 1 3 0", reordering = 3:1)
##!}

`inla.matrix2graph` = function(...)
{
    ## this function is not really needed, but is included for
    ## backwards compatibility
    return (inla.read.graph(...))
}

`inla.graph2matrix` = function(...)
{
    g = inla.read.graph(...)
    i = rep(1:g$n,  1L+g$nnbs)
    j = unlist(sapply(1:g$n, function(i,g) return (c(i, g$nbs[[i]])), g))
    stopifnot(length(i) == length(j))
    Q = inla.as.dgTMatrix(sparseMatrix(i = i,  j = j,  x = rep(1, length(i))))

    return (Q)
}

`inla.spy` = function(..., reordering = NULL)
{
    Q = inla.graph2matrix(...)
    Q[ Q != 0 ] = 1L
    Q  = 1L -Q
    if (!is.null(reordering)) {
        ## I don't know why this code fail *sometimes*:
        ##     Q[reordering, reordering] = Q
        ## giving this error:
        ##     Error in which(sel)[!vN0] : invalid subscript type 'S4'".
        ## So we need to work around the problem defining the inverse reordering instead
        r = reordering
        r[reordering] = 1:length(r)
        Q = Q[r, r]
    }

    inla.display.matrix(Q, nlevel = 2)
}
