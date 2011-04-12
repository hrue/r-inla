\name{inla.mesh.basis}
\alias{inla.mesh.basis}
\title{
  Basis functions for an inla.mesh
}
\description{
  Calculate basis functions on an \code{\link{inla.mesh}}
}
\usage{
inla.mesh.basis(mesh,
                type=c("b.spline", "sph.harm"),
                n=3,
                degree=2,
                knot.placement=c("uniform.area","uniform.latitude"),
                rot.inv=TRUE)
}
\arguments{
  \item{mesh}{
  }
  \item{type}{
  }
  \item{n}{
    For B-splines, the number of basis functions.  For spherical
    harmonics, \code{n} is the maximal harmonic order.
  }
  \item{degree}{
    Degree of B-spline polynomials.
  }
  \item{knot.placement}{
  }
  \item{rot.inv}{
  }
}
\details{
}
\value{
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
  \code{\link{inla.mesh.create}}
}
\examples{
n = 100
loc = matrix(runif(n*2), n, 2)
mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
basis = inla.mesh.basis(mesh, n=5)

proj = inla.mesh.projector(mesh)
image(proj$x, proj$y, inla.mesh.projector(proj, basis[,2]))
plot(old.mesh.class(mesh), basis[,2], draw.edges=FALSE, draw.vertices=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
