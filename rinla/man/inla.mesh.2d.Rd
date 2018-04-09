\name{inla.mesh.2d}
\alias{inla.mesh.2d}
\title{
High-quality triangulations
}
\description{
  Create a triangle mesh based on initial point locations, specified or
  automatic boundaries, and mesh quality parameters.
}
\usage{
inla.mesh.2d(loc = NULL,
             loc.domain = NULL,
             offset = NULL,
             n = NULL,
             boundary = NULL,
             interior = NULL,
             max.edge,
             min.angle = NULL,
             cutoff = 1e-12,
             max.n.strict = NULL,
             max.n = NULL,
             plot.delay = NULL,
             crs = NULL)
}
\arguments{
  \item{loc}{
    Matrix of point locations to be used as initial triangulation nodes.
    Can alternatively be a \code{SpatialPoints} or \code{SpatialPointsDataFrame} object.
  }
  \item{loc.domain}{
    Matrix of point locations used to determine the domain extent.
    Can alternatively be a \code{SpatialPoints} or \code{SpatialPointsDataFrame} object.
  }
  \item{offset}{
    The automatic extension distance.  One or two values, for an inner
    and an optional outer extension.  If negative, interpreted as a
    factor relative to the approximate data diameter (default=-0.10???)
  }
  \item{n}{
    The number of initial nodes in the automatic extensions (default=16)
  }
  \item{boundary}{
    A list of one or two \code{\link{inla.mesh.segment}} objects
    describing domain boundaries.
  }
  \item{interior}{
    An  \code{\link{inla.mesh.segment}} object describing desired
    interior edges.
  }
  \item{max.edge}{
    The largest allowed triangle edge length.  One or two values.
  }
  \item{min.angle}{
    The smallest allowed triangle angle.  One or two values. (Default=21)
  }
  \item{cutoff}{
    The minimum allowed distance between points.  Point at most as far
    apart as this are replaced by a single vertex prior to the mesh
    refinement step.
  }
  \item{max.n.strict}{The maximum number of vertices allowed,
    overriding \code{min.angle} and \code{max.edge}
    (default=-1, meaning no limit).
    One or two values, where the second value gives the number of
    additional vertices allowed for the extension.
  }
  \item{max.n}{The maximum number of vertices allowed,
    overriding \code{max.edge} only (default=-1, meaning no limit).
    One or two values, where the second value gives the number of
    additional vertices allowed for the extension.
  }
  \item{plot.delay}{
    On Linux (and Mac if appropriate X11 libraries are installed),
    specifying a nonnegative numeric value activates a rudimentary
    plotting system in the underlying \code{fmesher} program, showing
    the triangulation algorithm at work, with waiting time factor
    \code{plot.delay} between each step.

    On all systems, specifying any negative value activates displaying
    the result after each step of the multi-step domain extension
    algorithm.
  }
  \item{crs}{
    An optional \code{CRS} or \code{inla.CRS} object
  }
}
%\details{
  %%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An \code{inla.mesh} object.
%%  ~Describe the value returned
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
\code{\link{inla.mesh.create}},
\code{\link{inla.delaunay}},
\code{\link{inla.nonconvex.hull}}
}
\examples{
loc <- matrix(runif(10*2),10,2)

if (require("splancs")) {
  boundary <- list(inla.nonconvex.hull(loc, 0.1, 0.15),
                   inla.nonconvex.hull(loc, 0.2, 0.2))
  offset <- NULL
} else {
  boundary <- NULL
  offset <- c(0.1, 0.2)
}
mesh <- inla.mesh.2d(loc, boundary=boundary, offset=offset, max.edge=c(0.05, 0.1))

plot(mesh)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
