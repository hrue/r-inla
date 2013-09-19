\name{inla.mesh.lattice}
\alias{inla.mesh.lattice}
\title{
  Lattice grids for inla.mesh
}
\description{
  Construct a lattice grid for \code{\link{inla.mesh}}
}
\usage{
inla.mesh.lattice(x = seq(0, 1, length.out=2),
          y = seq(0, 1, length.out=2),
          z = NULL,
          dims = if (is.matrix(x)) {
                   dim(x)
                 } else {
                   c(length(x), length(y))
                 },
          units = NULL)
}
\arguments{
  \item{x}{
  }
  \item{y}{
  }
  \item{z}{
  }
  \item{dims}{
  }
  \item{units}{
    One of \code{c("default", "longlat", "longsinlat")}.
  }
}
%\details{
%}
\value{
  An \code{inla.mesh.lattice} object.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  \code{\link{inla.mesh}}
}
\examples{
lattice = inla.mesh.lattice(seq(0, 1, length.out=17), seq(0, 1, length.out=10))

## Use the lattice "as-is", without refinement:
mesh = inla.mesh.create(lattice=lattice, boundary=lattice$segm)
mesh = inla.mesh.create(lattice=lattice, extend=FALSE)
plot(mesh)

## Refine the triangulation, with limits on triangle angles and edges:
mesh = inla.mesh.create(lattice=lattice,
                        refine=list(max.edge=0.08),
                        extend=FALSE)
plot(mesh)

## Add an extension around the lattice, but maintain the lattice edges:
mesh = inla.mesh.create(lattice=lattice,
                        refine=list(max.edge=0.08),
                        interior=lattice$segm)
plot(mesh)

## Only add extension:
mesh = inla.mesh.create(lattice=lattice, refine=list(max.edge=0.08))
plot(mesh)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
