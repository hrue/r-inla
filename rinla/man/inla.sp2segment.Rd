\name{inla.sp2segment}
\alias{inla.sp2segment.Polygon}
\alias{inla.sp2segment.Polygons}
\alias{inla.sp2segment.SpatialPolygons}
\alias{inla.sp2segment.SpatialPolygonsDataFrame}
\title{
  Convert \code{sp} curve objects to \code{inla.mesh.segment} objects.
}
\description{
  Convert \code{sp} curve objects to \code{inla.mesh.segment} objects.
}
\usage{
inla.sp2segment(...)

\method{inla.sp2segment}{Polygon}(sp, ...)
\method{inla.sp2segment}{Polygons}(sp, join=TRUE, ...)
\method{inla.sp2segment}{SpatialPolygons}(sp, join=TRUE, grp=NULL, ...)
\method{inla.sp2segment}{SpatialPolygonsDataFrame}(sp, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{
    An \code{sp} polygon object of class \code{Polygon},
    \code{Polygons}, \code{SpatialPolygons}, or
    \code{SpatialPolygonsDataFrame}.
  }
  \item{join}{
    If \code{TRUE}, join multiple polygons into a single segment (possibly
    non-simply connected).
  }
  \item{grp}{
    Group ID specification for each polygon, as used by
    \code{\link{inla.mesh.segment}}, one ID per polygon.
}
\item{\dots}{
  Additional arguments passed on to other methods.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A \code{\link{inla.mesh.segment}} object, or a list of
  \code{\link{inla.mesh.segment}} objects.
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
\code{\link{inla.mesh.segment}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
