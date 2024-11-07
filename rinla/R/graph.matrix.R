#' Construct a neighbour-matrix from a `graph`
#' 
#' Construct a neighbour-matrix from a `graph` and disaply it
#' 
#' 
#' @aliases graph2matrix inla.graph2matrix spy inla.spy
#' @param graph An `inla.graph`-object, a (sparse) symmetric matrix, a
#' filename containing the graph, or a list or collection of characters and/or
#' numbers defining the graph.
#' @param reordering A possible reordering. Typical the one obtained from a
#' `inla`-call, `result$misc$reordering`, or the result of
#' `inla.qreordering`.
#' @param ... Additional arguments to `inla.read.graph()`
#' @return `inla.graph2matrix` returns a sparse symmetric matrix where the
#' non-zero pattern is defined by the `graph`.  The `inla.spy`
#' function, plots the associated sparse matrix defined by the `graph`. The `reordering`
#' argument is typically the reordering found by `inla.qreordering()`.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso [inla.read.graph()], [inla.qreordering()]
#' @examples
#' 
#' n = 50
#' Q = matrix(0, n, n)
#' idx = sample(1:n, 2*n, replace=TRUE)
#' Q[idx, idx] = 1
#' diag(Q) = 1
#' g = inla.read.graph(Q)
#' QQ = inla.graph2matrix(g)
#' inla.spy(QQ)
#' print(all.equal(as.matrix(Q), as.matrix(QQ)))
#' 
#' g.file = inla.write.graph(g, filename = tempfile())
#' inla.dev.new()
#' inla.spy(g.file)
#' inla.spy(g.file,  reordering = inla.qreordering(g))
#' 
#' g = inla.read.graph(g.file)
#' inla.dev.new()
#' inla.spy(g)
#' 
#' ## Old examples that don't work with the inla.spy call syntax:
#' # inla.dev.new()
#' # inla.spy(3, 1, "1 2 2 1 1 3 0")
#' # inla.dev.new()
#' # inla.spy(3, 1, "1 2 2 1 1 3 0", reordering = 3:1)
#' 
#' @name graph.matrix
#' @rdname graph.matrix
#' @export

`inla.matrix2graph` <- function(graph, ...) {
    ## this function is not really needed, but is included for
    ## backwards compatibility
    return(inla.read.graph(graph, ...))
}

#' @rdname graph.matrix
#' @export
`inla.graph2matrix` <- function(graph, ...) {
    g <- inla.read.graph(graph, ...)
    i <- rep(1:g$n, 1L + g$nnbs)
    j <- unlist(sapply(1:g$n, function(i, g) {
        return(c(i, g$nbs[[i]]))
    }, g))
    stopifnot(length(i) == length(j))
    Q <- inla.as.dgTMatrix(sparseMatrix(i = i, j = j, x = rep(1, length(i))))

    return(Q)
}

#' @rdname graph.matrix
#' @export
`inla.spy` <- function(graph, ..., reordering = NULL) {
    ## add this test here, as otherwise, this can be very inefficient
    ## for large matrices. this is because we convert it into a graph
    ## and then back to a matrix.
    M <- try(inla.as.sparse(graph), silent = TRUE)
    if (inherits(M, "try-error")) {
        M <- inla.graph2matrix(graph, ...)
    }
    if (!is.null(reordering)) {
        if (any(names(reordering) %in% "reordering")) {
            r <- reordering$reordering
        } else {
            r <- reordering
        }
        M@i <- as.integer(r[1L + M@i]) - 1L
        M@j <- as.integer(r[1L + M@j]) - 1L
    }
    d <- dim(M)
    xx <- c(0, d[2]+1, d[2]+1, 0, 0)
    yy <- c(0, 0, d[1]+1, d[1]+1, 0)
    plot(xx, yy, type = "l", bty = "L",
         xlim = c(0, d[2]+1), ylim = c(0, d[1]+1), lwd = 2, asp = 1)
    points(1L+M@j, d[1] - M@i, pch = 19, ...)
}
