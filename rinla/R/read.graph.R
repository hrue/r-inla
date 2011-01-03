## This function interface with mode = graph in the inlaprogram.

##!\name{read.graph}
##!\alias{read.graph}
##!\alias{inla.read.graph}
##!\title{Read a graph-file}
##!\description{Reads a graph spesification from file and create an R-object for it, and write a graph to file."}
##!\usage{
##!inla.read.graph(graph.file)
##!inla.write.graph(graph, graph.file)}
##!\arguments{
##!    \item{graph.file}{The filename of the graph.}
##!    \item{graph}{An graph-object (output from \code{inla.read.graph}).}
##!}
##!\value{
##!    The output of \code{inla.read.graph} is a graph-object, \code{graph} say, where
##!    \code{graph$n} is the size of the graph,
##!    \code{graph$nnbs} is a vector with the number of neigbours, 
##!    \code{graph$nbs} is a list-list with the neigbours, and
##!    \code{graph$graph.file} is the filename for which this graph is based on.
##!    \code{graph$cc$id} is a vector with the connected component id for each node (starting from 1)
##!    \code{graph$cc$n}  is the number of connected components
##!    \code{graph$cc$nodes}  is a list-list of nodes belonging to each connected component
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{inla.debug.graph}
##!\examples{
##!cat("2 1 1 2 2 1 1\n", file="g.dat")
##!g = inla.read.graph("g.dat")
##!print(g)
##!
##!## this function just writes the graph back to file; useful
##!## for converting between a 0-based and 1-based graph.
##!inla.write.graph(g, "gg.graph")
##!}

`inla.read.graph` = function(graph.file)
{
    if (missing(graph.file))
        stop("Need 'graph.file'")

    stopifnot(file.exists(graph.file))
              
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.call.builtin()), "-s -m graph", graph.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.call.builtin()), "-s -m graph", graph.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    ## since several nodes occure in the same line... we get output
    ## like "2" "0 1 1" "1 1 0" etc...
    s = as.integer(unlist(sapply(s, function(x) strsplit(x, " "))))

    n = s[1]
    g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), graph.file = graph.file,
            cc = list(id = NA, n = NA, nodes = NA))
    k = 2
    for(i in 1L:n) {
        
        stopifnot(s[k]+1L == i)
        k = k+1

        g$nnbs[i] = s[k]
        k = k+1

        if (g$nnbs[i] > 0) {
            g$nbs[[i]] = s[k:(k + g$nnbs[i] -1)] + 1L
            k = k + g$nnbs[i]
        }
    }
    stopifnot(k - 1L + n == length(s))

    ## admin the connected components-info
    g$cc$id = as.integer(s[k:length(s)]) + 1L # 0-based
    g$cc$n = max(g$cc$id)
    g$cc$nodes = lapply(1:g$cc$n, function(cc.id, id) which(cc.id == id), cc = g$cc$id)

    return (g)
}
`inla.write.graph` = function(graph, graph.file = "graph.dat") {
    ## write a graph read from inla.read.graph, or in that format, to
    ## file.
    fd = file(graph.file , "w")
    cat(graph$n, "\n", file = fd)
    for(i in 1:graph$n) {
        cat(i, graph$nnbs[i], graph$nbs[[i]], "\n", file = fd)
    }
    close(fd)

    return (invisible())
}

