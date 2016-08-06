\name{as.inla.mesh.segment}
\alias{as.inla.mesh.segment}
\alias{as.inla.mesh.segment.Line}
\alias{as.inla.mesh.segment.Lines}
\alias{as.inla.mesh.segment.SpatialLines}
\alias{as.inla.mesh.segment.SpatialLinesDataFrame}
\alias{as.inla.mesh.segment.Polygon}
\alias{as.inla.mesh.segment.Polygons}
\alias{as.inla.mesh.segment.SpatialPolygons}
\alias{as.inla.mesh.segment.SpatialPolygonsDataFrame}
\alias{inla.sp2segment}
\title{
  Convert \code{sp} curve objects to \code{inla.mesh.segment} objects.
}
\description{
  Convert \code{sp} curve objects to \code{inla.mesh.segment} objects.
}
\usage{
as.inla.mesh.segment(sp, ...)
inla.sp2segment(sp, ...) ## For backwards compatibility

\method{as.inla.mesh.segment}{Line}(sp, reverse=FALSE, crs=NULL, ...)
\method{as.inla.mesh.segment}{Lines}(sp, join=TRUE, crs=NULL, ...)
\method{as.inla.mesh.segment}{SpatialLines}(sp, join=TRUE, grp=NULL, ...)
\method{as.inla.mesh.segment}{SpatialLinesDataFrame}(sp, ...)
\method{as.inla.mesh.segment}{Polygon}(sp, crs=NULL, ...)
\method{as.inla.mesh.segment}{Polygons}(sp, join=TRUE, crs=NULL, ...)
\method{as.inla.mesh.segment}{SpatialPolygons}(sp, join=TRUE, grp=NULL, ...)
\method{as.inla.mesh.segment}{SpatialPolygonsDataFrame}(sp, ...)
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
  \item{reverse}{
    Logical, indicating if the line sequence should be traversed backwards.
  }
  \item{crs}{
    An optional \code{CRS} or \code{inla.CRS} object
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
