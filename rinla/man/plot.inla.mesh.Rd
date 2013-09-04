\name{plot.inla.mesh}
\alias{plot.inla.mesh}
\title{
  Draw a triangulation mesh object
}
\description{
  Plots an \code{\link{inla.mesh}} object using either standard graphics
  or with \code{rgl}.
}
\usage{
\method{plot}{inla.mesh}(x,
     col = "white",
     t.sub = 1:nrow(mesh$graph$tv),
     add = FALSE,
     lwd = 1,
     xlim = range(mesh$loc[, 1]),
     ylim = range(mesh$loc[, 2]),
     main = NULL,
     rgl = FALSE,
     size = 2,
     draw.vertices = FALSE,
     vertex.color = "black",
     draw.edges = TRUE,
     edge.color = rgb(0.3, 0.3, 0.3),
     draw.segments = draw.edges, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An \code{\link{inla.mesh}} object.
}
\item{col}{
  Color specification.  A single named color, a vector of scalar values,
  or a matrix of RGB values.  Requires \code{rgl=TRUE}.
}
\item{t.sub}{
  Optional triangle index subset to be drawn.
}
\item{add}{
  If \code{TRUE}, adds to the current plot instead of starting a new one.
}
\item{lwd}{
  Line width for triangle edges.
}
\item{xlim}{
  X-axis limits.
}
\item{ylim}{
  Y-axis limits.
}
\item{main}{
  The main plot title.  If not specified, a default title is generated
  based on the mesh type.
}
\item{rgl}{
  When \code{TRUE}, generates an \code{rgl} plot instead of a generic
  graphics plot.  Allows 3D plotting and color surface plotting.
}
\item{size}{
  Size of vertex points in \code{rgl} plotting. See \code{rgl.material}.
}
\item{draw.vertices}{
  If \code{TRUE}, draw triengle vertices.
}
  \item{vertex.color}{
  Color specification for all vertices.
}
  \item{draw.edges}{
  If \code{TRUE}, draw triangle edges.
}
\item{edge.color}{
  Color specification for all edges.
}
  \item{draw.segments}{
  If \code{TRUE}, draw boundary and interior constraint edges more
  prominently.
}
\item{\dots}{
  Further graphics parameters, interpreted by the respective plotting
  systems.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
\code{\link{plot.inla.trimesh}}
}
\examples{
mesh = inla.mesh.create(globe=10)
plot(mesh)
\donttest{
if (require(rgl)) {
  plot(mesh, rgl=TRUE, col=mesh$loc[,1])
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
