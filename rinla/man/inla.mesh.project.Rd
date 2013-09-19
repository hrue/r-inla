\name{inla.mesh.project}
\alias{inla.mesh.project}
\alias{inla.mesh.projector}
\alias{inla.mesh.projector.inla.mesh}
\alias{inla.mesh.project.inla.mesh}
\alias{inla.mesh.project.inla.mesh.projector}
\alias{inla.mesh.project.inla.mesh.1d}
\alias{inla.mesh.projector.inla.mesh.1d}
\title{
  Methods for projecting to/from an inla.mesh
}
\description{
  Calculate a lattice projection to/from an \code{\link{inla.mesh}}
}
\usage{
inla.mesh.project(...)
inla.mesh.projector(...)

\method{inla.mesh.projector}{inla.mesh}(mesh,
             loc = NULL,
             lattice = NULL,
             xlim = range(mesh$loc[,1]),
             ylim = range(mesh$loc[,2]),
             dims = c(100,100),
             projection = NULL,
             ...)

\method{inla.mesh.projector}{inla.mesh.1d}(mesh,
             loc = NULL,
             xlim = mesh$interval,
             dims = 100, ...)

\method{inla.mesh.project}{inla.mesh.projector}(projector, field, ...)

\method{inla.mesh.project}{inla.mesh}(mesh, loc, field = NULL, ...)
\method{inla.mesh.project}{inla.mesh.1d}(mesh, loc, field = NULL, ...)
}
\arguments{
  \item{mesh}{
    An \code{\link{inla.mesh}} or \code{\link{inla.mesh.1d}} object.
  }
  \item{loc}{
    Projection locations.
  }
  \item{lattice}{
    An \code{\link{inla.mesh.lattice}} object.
  }
  \item{xlim}{
    X-axis limits for a lattice.
  }
  \item{ylim}{
    Y-axis limits for a lattice.
  }
  \item{dims}{
    Lattice dimensions.
  }
  \item{projector}{An \code{inla.mesh.projector} object.}
  \item{field}{
    Basis function weights, one per mesh basis function, describing the
    function to be avaluated at the projection locationssFunction values for on the mesh 
  }
  \item{projection}{
    One of \code{c("default", "longlat", "longsinlat", "mollweide")}.
  }
  \item{\dots}{Additional arguments passed on to methods.}
}
\details{
  The call \code{inla.mesh.project(mesh, loc, field=..., ...)}, is a
  shortcut to inla.mesh.project(inla.mesh.projector(mesh, loc), field).
}
\value{
  For \code{inla.mesh.project(mesh, ...)}, a list with projection information.
  For \code{inla.mesh.projector(mesh, ...)}, an \code{inla.mesh.projector} object.
  For \code{inla.mesh.project(projector, field, ...)}, a field projected
  from the mesh onto the locations given by the projector object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh}},
  \code{\link{inla.mesh.1d}},
  \code{\link{inla.mesh.lattice}}
}
\examples{
n = 20
loc = matrix(runif(n*2), n, 2)
mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
proj = inla.mesh.projector(mesh)
field = cos(mesh$loc[,1]*2*pi*3)*sin(mesh$loc[,2]*2*pi*7)
image(proj$x, proj$y, inla.mesh.project(proj, field))
\donttest{
if (require(rgl)) {
  plot(mesh, rgl=TRUE, col=field, draw.edges=FALSE, draw.vertices=FALSE)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
