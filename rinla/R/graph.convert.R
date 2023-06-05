## various graph convert utilities that I got from someone. I don't
## know how useful they are really...

#' INLA utility functions
#' 
#' Various utility functions for INLA
#' 
#' 
#' @aliases geobugs2inla inla.geobugs2inla
#' @param adj A vector listing the ID numbers of the adjacent areas for each
#' area. This is a sparse representation of the full adjacency matrix for the
#' study region, and can be generated using the Adjacency Tool from the Map
#' menu in GeoBUGS.
#' @param num A vector of length N (the total number of areas) giving the
#' number of neighbours n.i for each area.
#' @param graph.file Name of the file of the new graph in the INLA format.
#' @return The return value is the name of the graph-file created.
#' @note These are all the same function, and the two different names are due
#' to backward-compatibility
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso [inla()], [inla.surv()],
#' [hyperpar.inla()]
#' @keywords graph
#' 
#' @name graph.convert
#' @rdname graph.convert
#' @export

`inla.geobugs2inla` <- function(adj, num, graph.file = "graph.dat") {
    return(inla.graph.convert.2(adj, num, graph.file))
}

`inla.graph.convert.2` <- function(adj, num, graph.file = "graph.dat") {
    ## A function for converting GeoBUGS adjacency data into the INLA
    ## graph format. Kindly provided by Aki Havunlinna tkk.fi; thanks.

    fd <- file(graph.file, "w")
    len <- length(num)
    cat(len, "\n", file = fd)
    k <- 1L
    for (i in 1L:len) {
        if (num[i] > 0) {
            cat(i, num[i], adj[k:(k + num[i] - 1L)], "\n", file = fd)
        } else {
            cat(i, num[i], "\n", file = fd)
        }
        k <- k + num[i]
    }
    close(fd)
    return(graph.file)
}
