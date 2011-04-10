\name{read.graph}
\alias{read.graph}
\alias{write.graph}
\alias{inla.read.graph}
\alias{inla.write.graph}
\title{Read a graph-file}
\description{Reads a graph spesification from file and create an R-object for it, and write a graph to file.}
\usage{
graph = inla.read.graph(graph.file)
inla.write.graph(graph, graph.file)}
\arguments{
    \item{graph.file}{The filename of the graph.}
    \item{graph}{An graph-object (output from \code{inla.read.graph}).}
}
\value{
    The output of \code{inla.read.graph} is a graph-object, \code{graph} say, where
    \code{graph$n} is the size of the graph,
    \code{graph$nnbs} is a vector with the number of neigbours, 
    \code{graph$nbs} is a list-list with the neigbours, and
    \code{graph$graph.file} is the filename for which this graph is based on.
    \code{graph$cc$id} is a vector with the connected component id for each node (starting from 1)
    \code{graph$cc$n}  is the number of connected components
    \code{graph$cc$nodes}  is a list-list of nodes belonging to each connected component

    Methods implemented for \code{inla.graph} is \code{summary}, \code{print} and \code{plot}.
}
\author{Havard Rue \email{hrue@math.ntnu.no}}
\seealso{inla.debug.graph}
\examples{
cat("2 1 1 2 2 1 1\n", file="g.dat")
g = inla.read.graph("g.dat")
print(g)
plot(g)

## this function just writes the graph back to file; useful
## for converting between a 0-based and 1-based graph.
inla.write.graph(g, "gg.graph")
}
