\name{mesh.segm}
\alias{mesh.segm}
\alias{inla.mesh.segm}
\title{
  Constraint segments for inla.mesh
}
\description{
  Constructs an \code{inla.mesh.segm} object that can be used to specify
  boundary and interior constraint edges in calls to \code{\link{inla.mesh}}.
}
\usage{
mesh.segm(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE)
}
\arguments{
  \item{loc}{
    Matrix of point locations.
  }
  \item{idx}{
    Segment index sequence vector or index pair matrix.  The indices
    refer to the rows of \code{loc}.  If \code{loc==NULL}, the indices
    will be interpreted as indices into the point specification supplied
    to \code{\link{inla.mesh}}.
    If \code{is.bnd==TRUE}, defaults to linking all the points in
    \code{loc}, as \code{c(1:nrow(loc),1L)}, otherwise \code{1:nrow(loc)}.  
  }
  \item{grp}{
    Vector of group labels for each segment.
    Set to \code{NULL} to let the labels be chosen automatically in a
    call to \code{\link{inla.mesh}}.
  }
  \item{is.bnd}{
    \code{TRUE} if the segments are boundary segments, otherwise \code{FALSE}.
  }
}
%%\details{
%%  Some text. \describe{Some more text} Further text.
%%}
\value{
  An \code{inla.mesh.segm} object.
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
segm.bnd = mesh.segm(loc.bnd)
segm.int = mesh.segm(loc.int, is.bnd=FALSE)

## Points to be meshed
loc = matrix(runif(10*2),10,2)*0.9+0.05
mesh = inla.mesh(loc, boundary=segm.bnd, interior=segm.int,
                 refine=list())
plot(mesh)

\dontrun{
mesh = inla.mesh(loc, interior=list(segm.bnd, segm.int))
plot(mesh)
}
}
%\keyword{ fmesher }
