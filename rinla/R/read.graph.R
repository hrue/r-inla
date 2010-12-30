## This function interface with mode = graph in the inlaprogram.

##!\name{read.graph}
##!\alias{read.graph}
##!\alias{inla.read.graph}
##!\title{Read a graph-file}
##!\description{Reads a graph spesification from file and create an R-object for it"}
##!\usage{inla.read.graph(graph.file)}
##!\arguments{
##!    \item{graph.file}{The filename of the graph.}
##!}
##!\value{
##!    The output is a graph-object, \code{graph} say, where
##!    \code{graph$n} is the size of the graph,
##!    \code{graph$nnbs} is a vector with the number of neigbours, 
##!    \code{graph$nbs} is a list-list with the neigbours, and
##!    \code{graph$graph.file} is the filename for which this graph is based on.
##!    \code{graph$cc$id} is a vector with the connected component id for each node (starting from 1)
##!    \code{graph$cc$n}  is the number of connected components
##!    \code{graph$cc$nodes}  is a list-list of nodes belonging to each connected component
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{}
##!\examples{
##!cat("2 1 1 2 2 1 1\n", file="g.dat")
##!g = inla.read.graph("g.dat")
##!print(g)
##!}

`inla.read.graph` = function(graph.file)
{
    if (missing(graph.file))
        stop("Need 'graph.file'")

    stopifnot(file.exists(graph.file))
              
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m graph", graph.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m graph", graph.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    ## since several nodes occure in the same line... we get output
    ## like "2" "0 1 1" "1 1 0" etc...
    s = as.integer(unlist(sapply(s, function(x) strsplit(x, " "))))

    n = s[1]
    g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), graph.file = graph.file, cc = list(id = NA, n = NA, nodes = NA))

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
