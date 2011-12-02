## WARNING! DO NOT EDIT!
## File automatically generated from read.graph.R
##
\name{read.graph}
\alias{read.graph}
\alias{write.graph}
\alias{inla.read.graph}
\alias{inla.write.graph}
\alias{inla.graph}
\title{Read a graph-file}
\description{Reads a graph specification from file and create an R-object for it, and write a graph to file.}
\usage{
graph = inla.read.graph(graph.file)
inla.write.graph(graph, graph.file)}
\arguments{
    \item{graph.file}{The filename of the graph.}
    \item{graph}{An \code{inla.graph} object (output from \code{inla.read.graph}).}
}
\value{
    The output of \code{inla.read.graph} is an \code{inla.graph} object, with elements
    \item{n }{is the size of the graph}
    \item{nnbs }{is a vector with the number of neigbours}
    \item{nbs }{is a list-list with the neigbours}
    \item{graph.file }{is the name of the file on which this graph is based}
    \item{cc }{list with connected component information
        \itemize{
            \item{\code{id} }{is a vector with the connected component id for each node (starting from 1)}
            \item{\code{n} }{is the number of connected components}
            \item{\code{nodes} }{is a list-list of nodes belonging to each connected component}
        }
    }
    Methods implemented for \code{inla.graph} are \code{summary}, \code{print} and \code{plot}.
}
\author{Havard Rue \email{hrue@math.ntnu.no}}
\seealso{
    \code{\link{inla.debug.graph}},
    \code{\link{inla.matrix2graph}}
}
\examples{
cat("2 1 1 2 2 1 1\n", file="g.dat")
g = inla.read.graph("g.dat")
print(g)
plot(g)

## this function just writes the graph back to file; useful
## for converting between a 0-based and 1-based graph.
inla.write.graph(g, "gg.graph")
}
