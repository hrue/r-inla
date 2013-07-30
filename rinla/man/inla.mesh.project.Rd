\name{inla.mesh.project}
\alias{inla.mesh.project}
\alias{inla.mesh.projector}
\alias{inla.mesh.projector.inla.mesh}
\alias{inla.mesh.project.inla.mesh}
\alias{inla.mesh.project.inla.mesh.projector}
\title{
  Methods for projecting to/from an inla.mesh
}
\description{
  Calculate a lattice projection to/from an \code{\link{inla.mesh}}
}
\usage{
inla.mesh.project(mesh, loc)
inla.mesh.projector(mesh, xlim, ylim, dims, projection)
inla.mesh.project(projector, field)
}
\arguments{
  \item{mesh}{
  }
  \item{loc}{
  }
  \item{xlim}{
  }
  \item{ylim}{
  }
  \item{dims}{
  }
  \item{projection}{
    One of \code{c("default", "longlat", "longsinlat")}.
  }
}
\details{
}
\value{
  For \code{inla.mesh.project(mesh, ...)}, a list with projection information.
  For \code{inla.mesh.projector(mesh, ...)}, an \code{inla.mesh.projector} object.
  For \code{inla.mesh.project(projector, ...)}, a field projected
  from the mesh onto a grid.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh}}
}
\examples{
n = 20
loc = matrix(runif(n*2), n, 2)
mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
proj = inla.mesh.projector(mesh)
field = cos(mesh$loc[,1]*2*pi*3)*sin(mesh$loc[,2]*2*pi*7)
image(proj$x, proj$y, inla.mesh.project(proj, field))
if (require(rgl)) {
  plot(mesh, rgl=TRUE, col=field, draw.edges=FALSE, draw.vertices=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
