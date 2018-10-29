## Export: inla.geobugs2inla

##!\name{geobugs2inla}
##!\alias{geobugs2inla}
##!\alias{inla.geobugs2inla}
##!\keyword{graph}
##!
##!\title{INLA utility functions}
##!
##!\description{Various utility functions for INLA}
##!
##!\usage{
##!     inla.geobugs2inla(adj, num, graph.file="graph.dat")
##!}
##!
##!\arguments{
##!
##!  \item{adj}{A vector listing the ID numbers of the adjacent areas for
##!  each area. This is a sparse representation of the full adjacency
##!  matrix for the study region, and can be generated using the Adjacency
##!  Tool from the Map menu in GeoBUGS.}
##!
##!  \item{num}{A vector of length N (the total number of areas) giving the
##!  number of neighbours n.i for each area.}
##!
##!  \item{graph.file}{Name of the file of the new graph in the INLA
##!    format.}
##!
##!}
##!
##!
##!\value{The return value is the name of the graph-file created.}
##!
##!\note{These are all the same function,  and the two different names are due to backward-compatibility}
##!
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!
##!\seealso{\code{\link{inla}}, \code{\link{inla.surv}}, \code{\link{hyperpar.inla}}}
##!
##!%%\examples{}

## various graph convert utilities that I got from someone. I don't
## know how useful they are really...

`inla.geobugs2inla` = function(adj, num, graph.file="graph.dat")
{
    return (inla.graph.convert.2(adj, num, graph.file))
}

`inla.graph.convert.2` = function(adj, num, graph.file="graph.dat")
{
    ## A function for converting GeoBUGS adjacency data into the INLA
    ## graph format. Kindly provided by Aki Havunlinna tkk.fi; thanks.

    fd = file(graph.file,  "w")
    len <- length(num)
    cat(len, '\n', file=fd)
    k = 1L
    for(i in 1L:len) {
        if (num[i] > 0) {
            cat(i, num[i], adj[k:(k+num[i]-1L)], "\n", file = fd)
        } else {
            cat(i, num[i], "\n", file=fd)
        }
        k = k + num[i]
    }
    close(fd)
    return (graph.file)
}
