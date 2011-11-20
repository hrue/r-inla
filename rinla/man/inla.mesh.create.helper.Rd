\name{inla.mesh.create.helper}
\alias{inla.mesh.create.helper}
\title{
  Helper for high-quality triangulations 
}
\description{
  Create a constrained refined Delaunay triangulation (CRDT) for a
  spatial domain.  Allows an interior domain with small triangles and an
  outer domain with large triangles.
}
\usage{
inla.mesh.create.helper(points = NULL,
                        points.domain = NULL,
                        offset = c(-0.05, -0.15),
                        n = c(8, 16),
                        boundary = NULL,
                        interior = NULL,
                        max.edge,
                        min.angle = c(21, 21),
                        cutoff = 0,
                        plot.delay = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{
%%     ~~Describe \code{points} here~~
}
  \item{points.domain}{
%%     ~~Describe \code{points.domain} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{boundary}{
%%     ~~Describe \code{boundary} here~~
}
  \item{interior}{
%%     ~~Describe \code{interior} here~~
}
  \item{max.edge}{
%%     ~~Describe \code{max.edge} here~~
}
  \item{min.angle}{
%%     ~~Describe \code{min.angle} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{plot.delay}{
%%     ~~Describe \code{plot.delay} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
  \code{\link{inla.spde.create}},
  \code{\link{inla.spde.make.A}}
}
\examples{
loc = matrix(runif(100*2)*1000,100,2)
mesh = inla.mesh.create.helper(points.domain=loc, max.edge=c(50,500))
boundary = inla.mesh.boundary(mesh)
interior = inla.mesh.interior(mesh)

loc = matrix(runif(10000*2)*1000,10000,2)
mesh = inla.mesh.create.helper(points=loc,
                               cutoff=50,
                               max.edge=c(50,500))
A = inla.spde.make.A(mesh, loc=loc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
