\name{inla.mesh.segment}
\alias{inla.mesh.segment}
\title{
  Constraint segments for inla.mesh
}
\description{
  Constructs an \code{inla.mesh.segment} object that can be used to specify
  boundary and interior constraint edges in calls to \code{\link{inla.mesh}}.
}
\usage{
inla.mesh.segment(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE)
}
\arguments{
  \item{loc}{
    Matrix of point locations.
  }
  \item{idx}{
    Segment index sequence vector or index pair matrix.  The indices
    refer to the rows of \code{loc}.  If \code{loc==NULL}, the indices
    will be interpreted as indices into the point specification supplied
    to \code{\link{inla.mesh.create}}.
    If \code{is.bnd==TRUE}, defaults to linking all the points in
    \code{loc}, as \code{c(1:nrow(loc),1L)}, otherwise \code{1:nrow(loc)}.  
  }
  \item{grp}{
    Vector of group labels for each segment.
    Set to \code{NULL} to let the labels be chosen automatically in a
    call to \code{\link{inla.mesh.create}}.
  }
  \item{is.bnd}{
    \code{TRUE} if the segments are boundary segments, otherwise \code{FALSE}.
  }
}
%%\details{
%%  Some text. \describe{Some more text} Further text.
%%}
\value{
  An \code{inla.mesh.segment} object.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
}
\references{
}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh}}
}
\examples{
## Create a square boundary and a diagonal interior segment
loc.bnd = matrix(c(0,0, 1,0, 1,1, 0,1), 4, 2, byrow=TRUE)
loc.int = matrix(c(0.9,0.1, 0.1,0.6), 2, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd)
segm.int = inla.mesh.segment(loc.int, is.bnd=FALSE)

## Points to be meshed
loc = matrix(runif(10*2),10,2)*0.9+0.05
mesh = inla.mesh.create(loc,
                        boundary=segm.bnd,
                        interior=segm.int,
                        refine=list())
plot(mesh)

\dontrun{
mesh = inla.mesh.create(loc, interior=list(segm.bnd, segm.int))
plot(mesh)
}
}
%\keyword{ fmesher }
