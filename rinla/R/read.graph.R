## This function interface with mode = graph in the inlaprogram.

##!\name{read.graph}
##!\alias{read.graph}
##!\alias{write.graph}
##!\alias{inla.read.graph}
##!\alias{inla.write.graph}
##!\alias{inla.graph}
##!\title{Read a graph-file}
##!\description{Reads a graph specification from file and create an R-object for it, and write a graph to file.}
##!\usage{
##!graph = inla.read.graph(graph.file)
##!inla.write.graph(graph, graph.file)}
##!\arguments{
##!    \item{graph.file}{The filename of the graph.}
##!    \item{graph}{An \code{inla.graph} object (output from \code{inla.read.graph}).}
##!}
##!\value{
##!    The output of \code{inla.read.graph} is an \code{inla.graph} object, with elements
##!    \item{n }{is the size of the graph}
##!    \item{nnbs }{is a vector with the number of neigbours}
##!    \item{nbs }{is a list-list with the neigbours}
##!    \item{graph.file }{is the name of the file on which this graph is based}
##!    \item{cc }{list with connected component information
##!        \itemize{
##!            \item{\code{id} }{is a vector with the connected component id for each node (starting from 1)}
##!            \item{\code{n} }{is the number of connected components}
##!            \item{\code{nodes} }{is a list-list of nodes belonging to each connected component}
##!        }
##!    }
##!    Methods implemented for \code{inla.graph} are \code{summary}, \code{print} and \code{plot}.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{
##!    \code{\link{inla.debug.graph}},
##!    \code{\link{inla.matrix2graph}}
##!}
##!\examples{
##!cat("2 1 1 2 2 1 1\n", file="g.dat")
##!g = inla.read.graph("g.dat")
##!print(g)
##!plot(g)
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
        s = system(paste(shQuote(inla.call.builtin()), "-s -m graph", shQuote(graph.file)), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.call.builtin()), "-s -m graph", shQuote(graph.file)), intern=TRUE)
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

    class(g) = "inla.graph"
    return (g)
}

`inla.write.graph` = function(graph, graph.file = "graph.txt")
{
    ## write a graph read from inla.read.graph, or in that format, to
    ## file.
    fd = file(graph.file , "w")
    cat(graph$n, "\n", file = fd)
    for(i in 1:graph$n) {
        cat(i, graph$nnbs[i], graph$nbs[[i]], "\n", file = fd)
    }
    close(fd)
    return (graph.file)
}

## I add here some tools to view and summarize a such graphs...
plot.inla.graph = function(
        graph,
        filter = c("neato", "fdp"),
        attrs = NULL,
        scale = 0.5,
        node.names = NULL,
        ...)
{
    require(Rgraphviz) || stop("Need library 'Rgraphviz' from the Bioconductor package")
    require(graph) || stop("Need library 'graph'")

    filter = match.arg(filter)
    if (is.null(attrs)) {
        attrs = getDefaultAttrs(layoutType = filter)
    }
    if (!is.null(node.names)) {
        stopifnot(length(node.names) == graph$n)
    } else {
        node.names = as.character(1:graph$n)
    }
    g <- new("graphNEL", nodes = node.names, edgemode = "undirected")
    for (i in 1L:graph$n) {
        if (graph$nnbs[i] > 0) {
            j = graph$nbs[[i]]
            j = j[j > i]
            if (length(j) > 0) {
                g = addEdge(node.names[i], node.names[j], g)
            }
        }
    }
    attrs$node$height = as.numeric(attrs$node$height) * scale
    attrs$node$width = as.numeric(attrs$node$width) * scale
    plot(g, filter, attrs = attrs, ...)
}

`summary.inla.graph` = function(graph, ...)
{
    ret = list()
    ret = c(ret, list(n = graph$n))
    ret = c(ret, list(ncc = graph$cc$n))
    ret = c(ret, list(nbs = table(graph$nnbs)))

    class(ret) = "inla.graph.summary"
    return(ret)
}

`print.inla.graph.summary` = function(go, ...)
{
    cat(paste("\tn = ",  go$n, "\n"))
    cat(paste("\tncc = ",  go$ncc, "\n"))
    cat(inla.paste(c("\tnbs = (names) ",  names(go$nbs), "\n")))
    cat(inla.paste(c("\t      (count) ",  go$nbs, "\n")))

    return(invisible())
}
