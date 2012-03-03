##!\name{graph2matrix}
##!\alias{graph2matrix}
##!\alias{inla.graph2matrix}
##!\alias{inla.matrix2graph}
##!\alias{spy}
##!\alias{inla.spy}
##!\title{Convert between a (possible) sparse symmetric matrix and an \code{inla.graph}-object}
##!\description{Convert between a (possible) sparse symmetric matrix and an \code{inla.graph}-object}
##!\usage{
##!Q = inla.graph2matrix(graph)
##!g = inla.matrix2graph(Q)
##!inla.spy(graph, ...)
##!\arguments{
##!    \item{graph}{An \code{inla.graph}-object, a filename containing an \code{inla.graph}-object,  or a (sparse) symmetric matrix}
##!    \item{g}{An \code{inla.graph}-object}
##!    \item{Q}{An (possible) sparse symmtric matrix}
##!}
##!\value{
##!     \code{inla.graph2matrix} returns a sparse symmetric matrix where the non-zero pattern is defined by the \code{graph}.
##!     The \code{inla.spy} function, plots a binary image of a \code{graph}.
##!
##!     The function\code{inla.read.graph}-function is the same as the function \code{inla.matrix2graph}.
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
##!inla.dev.new()
##!inla.spy(Q)
##!
##!g = inla.matrix2graph(Q)
##!inla.dev.new()
##!inla.spy(g)
##!
##!g.file = inla.write.graph(g)
##!inla.dev.new()
##!inla.spy(g.file)
##!}

`inla.matrix2graph` = function(...)
{
    return (inla.read.graph(...))
}

`inla.graph2matrix` = function(graph)
{
    if (missing(graph)) {
        return (NULL)
    }

    g = inla.read.graph(graph)
    i = rep(1:g$n,  1L+g$nnbs)
    j = unlist(sapply(1:g$n, function(i,g) return (c(i, g$nbs[[i]])), g))
    stopifnot(length(i) == length(j))
    Q = inla.as.dgTMatrix(sparseMatrix(i = i,  j = j,  x = rep(1, length(i))))

    return (Q)
}

`inla.spy` = function(graph, reordering = NULL, ...)
{
    if (missing(graph)) {
        return (NULL)
    }

    Q = inla.graph2matrix(graph)
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

    inla.display.matrix(Q, nlevel = 2, ...)
}
