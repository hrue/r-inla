\name{mesh.segm}
\alias{mesh.segm}
\title{
  Create a constraint segment object for fmesher
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mesh.segm(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE)
}
\arguments{
  \item{loc}{
    Matrix of point locations.
  }
  \item{idx}{
    Segment index sequence vector or index pair matrix.
    If \code{NULL}, defaults to linking all the points in \code{loc}, as
    \code{c(1:nrow(loc),1L)} (if \code{is.bnd==TRUE}) or
    \code{1:nrow(loc)} (if \code{is.bnd==FALSE})
  }
  \item{grp}{
    Vector of group labels for each segment.
    Set to \code{NULL} to let the labels to be chosen automatically in a
    call to \code{\link{inla.mesh}}.
  }
  \item{is.bnd}{
    \code{TRUE} if the segments are boundary segments, otherwise \code{FALSE}.
  }
}
\details{
  Some text. \describe{Some more text} Further text.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
Returns an \code{inla.fmesher.mesh.segm} object that can be used in calls
to \code{\link{inla.mesh}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create a square boundary and a diagonal interior segment
loc.bnd = matrix(c(0,0, 1,0, 1,1, 0,1), 4, 2, byrow=TRUE)
loc.int = matrix(c(0.9,0.1, 0.1,0.6), 2, 2, byrow=TRUE)
segm.bnd = mesh.segm(loc.bnd)
segm.int = mesh.segm(loc.int, is.bnd=FALSE)

## Points to be meshed
loc = matrix(runif(10*2),10,2)*0.9+0.05
mesh = inla.mesh(loc, interior=list(segm.bnd, segm.int))
mesh = inla.mesh(loc, boundary=segm.bnd, interior=segm.int,
                 refine=list())
}
\keyword{ fmesher }
