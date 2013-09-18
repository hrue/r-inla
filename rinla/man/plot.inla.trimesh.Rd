\name{plot.inla.trimesh}
\alias{plot.inla.trimesh}
\title{
  Low level triangulation mesh plotting
}
\description{
  Plots a triangulation mesh using \code{rgl}.
}
\usage{
\method{plot}{inla.trimesh}(x, S,
     color = NULL, color.axis = NULL,
     color.n = 512, color.palette = cm.colors, color.truncate = FALSE,
     alpha = NULL, lwd = 1, specular = "black",
     draw.vertices = TRUE,
     draw.edges = TRUE, edge.color = rgb(0.3, 0.3, 0.3),
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A 3-column triangle-to-vertex index map matrix.
}
\item{S}{
  A 3-column vertex coordinate matrix.
}
  \item{color}{
  Color specification.  A single named color, a vector of scalar values,
  or a matrix of RGB values.
}
  \item{color.axis}{
  The min/max limit values for the color mapping.
}
  \item{color.n}{
 The number of colors to use in the color palette.
}
  \item{color.palette}{
 A color palette function.
}
  \item{color.truncate}{
 If \code{TRUE}, truncate the colors at the color axis limits.
}
  \item{alpha}{
  Transparency/opaqueness values. See \code{rgl.material}.
}
\item{lwd}{
  Line width for edges. See \code{rgl.material}.
}
\item{specular}{
  Specular color. See \code{rgl.material}.
}
\item{draw.vertices}{
  If \code{TRUE}, draw triangle vertices.
}
  \item{draw.edges}{
  If \code{TRUE}, draw triangle edges.
}
\item{edge.color}{
  Edge color specification.
}
\item{\dots}{
  Additional parameters passed to and from other methods.
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
\code{\link{plot.inla.mesh}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
