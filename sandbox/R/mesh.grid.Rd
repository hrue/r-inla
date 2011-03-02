\name{mesh.grid}
\alias{mesh.grid}
\alias{inla.mesh.grid}
\title{
  Lattice grids for inla.mesh
}
\description{
  Construct a lattice grid for \code{\link{inla.mesh}}
}
\usage{
mesh.grid(x = seq(0, 1, length.out=2),
          y = seq(0, 1, length.out=2),
          z = NULL,
          dims = (inla.ifelse(is.matrix(x),
                              dim(x),
                              c(length(x), length(y)))),
          units = "default")
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
    One of \code{c("default", "longlat")}.
  }
}
\details{
}
\value{
  An \code{inla.mesh.grid} object.
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
grid = mesh.grid(seq(0, 1, length.out=17), seq(0, 1, length.out=10))

## Use the grid "as-is", without refinement:
mesh = inla.mesh(grid = grid, boundary=grid$segm)
mesh = inla.mesh(grid = grid, extend=FALSE)
plot(mesh)

\dontrun{
## Refine the triangulation, with limits on triangle angles and edges:
mesh = inla.mesh(grid = grid, refine=list(max.edge=0.08), extend=FALSE)
plot(mesh)

## Add an extension around the grid, but maintain the grid edges:
mesh = inla.mesh(grid = grid, refine=list(max.edge=0.08), interior=grid$segm)
plot(mesh)

## Only add extension:
mesh = inla.mesh(grid = grid, refine=list(max.edge=0.08))
plot(mesh)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
