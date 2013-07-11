\name{geobugs2inla}
\alias{geobugs2inla}
\alias{inla.util}
\alias{inla.utility}
\keyword{graph}

\title{INLA utility functions}

\description{Various utility functions for INLA}

\usage{
     geobugs2inla(adj, num, graph.file="graph.txt")
}

\arguments{

  \item{adj}{A vector listing the ID numbers of the adjacent areas for
  each area. This is a sparse representation of the full adjacency
  matrix for the study region, and can be generated using the Adjacency
  Tool from the Map menu in GeoBUGS.}

  \item{num}{A vector of length N (the total number of areas) giving the
  number of neighbours n.i for each area.}

  \item{graph.file}{Name of the file of the new graph in the INLA
    format.}

}


%%\value{}

%%\note{}

\author{Havard Rue \email{hrue@math.ntnu.no}}

\seealso{\code{\link{inla}}, \code{\link{inla.surv}}, \code{\link{hyperpar.inla}}}

%%\examples{}

