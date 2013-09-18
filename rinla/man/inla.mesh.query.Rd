\name{inla.mesh.query}
\alias{inla.mesh.query}
\title{
  High-quality triangulations
}
\description{
  Query information about an inla.mesh object.
}
\usage{
inla.mesh.query(mesh, ...)
}
\arguments{
  \item{mesh}{
    An \code{inla.mesh} object.
  }
  \item{\dots}{Query arguments.
    \itemize{
      \item{tt.neighbours}{
	Compute neighbour triangles for triangles; list of vectors:
	list(triangles, orders)
      }
      \item{vt.neighbours}{
	Compute neighbour triangles for vertices; list of vectors:
	list(vertices, orders)
      }
    }
  }
}
%\details{
%}
\value{
  A list of query results.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh.create}}, 
  \code{\link{inla.mesh.segment}}, \code{\link{inla.mesh.lattice}}
}
\examples{
loc = matrix(c(0.1,0.15),1,2)
lattice = inla.mesh.lattice(dims=c(10,10))
mesh = inla.mesh.create(loc=loc, lattice=lattice, extend=FALSE)

vt = which(inla.mesh.query(mesh,
                           vt.neighbours=list(mesh$idx$loc,
                                              4:6))$vt.neighbours)

mesh2 = inla.mesh.create(mesh$loc, tv=mesh$graph$tv[vt,,drop=FALSE],
                         refine=FALSE, extend=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
