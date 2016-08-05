\name{inla.mesh.segment}
\alias{inla.mesh.segment}
\alias{inla.mesh.segment.default}
\alias{inla.mesh.segment.inla.mesh.segment}
\alias{inla.contour.segment}
\title{
  Constraint segments for inla.mesh
}
\description{
  Constructs \code{inla.mesh.segment} objects that can be used to specify
  boundary and interior constraint edges in calls to
  \code{\link{inla.mesh}}.
}
\usage{
## Create or join inla.mesh.segment objects.
inla.mesh.segment(...)
\method{inla.mesh.segment}{default}(loc = NULL, idx = NULL, grp = NULL,
  is.bnd = TRUE, crs = NULL, ...)
\method{inla.mesh.segment}{inla.mesh.segment}(..., grp.default = 0)

inla.contour.segment(x = seq(0, 1, length.out = nrow(z)),
                     y = seq(0, 1, length.out = ncol(z)),
                     z,
                     nlevels = 10,
                     levels = pretty(range(z, na.rm = TRUE), nlevels),
                     groups = seq_len(length(levels)),
                     positive = TRUE,
                     eps = NULL,
                     crs = NULL)

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
  \item{grp.default}{
    When joining segments, use this group label for segments that have
    \code{grp=NULL}.
  }
  \item{x, y, z, nlevels, levels}{Parameters specifying a set of surface
    contours, with syntax described in \code{\link{contour}}.
}
  \item{groups}{Vector of group ID:s, one for each contour level.
}
  \item{positive}{\code{TRUE} if the contours should encircle positive
    level excursions in a counter clockwise direction.
}
\item{eps}{
  Tolerance for \code{\link{inla.simplify.curve}}.
}
\item{crs}{
  An optional \code{CRS} or \code{inla.CRS} object
}
\item{...}{
  Additional parameters.  When joining segments, a list of
  \code{inla.mesh.segment} objects.
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
%\references{
%}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh.create}},
  \code{\link{inla.mesh.2d}}
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
