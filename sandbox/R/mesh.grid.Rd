\name{mesh.grid}
\alias{mesh.grid}
\title{
  Construct a lattice grid for fmesher
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mesh.grid(x = seq(0, 1, length.out=2),
          y = seq(0, 1, length.out=2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An \code{fmesher.grid} object.
%%  ~Describe the value returned
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

grid = mesh.grid(seq(0, 1, length.out=17), seq(0, 1, length.out=10))

## Use the grid "as-is", without refinement:
mesh = inla.mesh(grid = grid, boundary=grid$segm)
mesh = inla.mesh(grid = grid, extend=FALSE)
plot(mesh)

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fmesher }
