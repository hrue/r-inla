\name{inla.nonconvex.hull}
\alias{inla.nonconvex.hull}
\alias{inla.nonconvex.hull.basic}
\title{
  Nonconvex set extensions.
}
\description{
  Constructs a nonconvex boundary for a point set using morphological
  operations.
}
\usage{
inla.nonconvex.hull(points,
                    convex = -0.15,
                    concave = convex,
                    resolution = 40,
                    eps = NULL)

inla.nonconvex.hull.basic(points,
                          convex = -0.15,
                          resolution = 40,
                          eps = NULL)
}
\arguments{
  \item{points}{
    2D point coordinates (2-column matrix).
}
\item{convex}{
  The desired extension radius.  Also determines the smallest allowed
  convex curvature radius.  Negative values are interpreted as fractions
  of the approximate initial set diameter.
}
\item{concave}{
  The desired minimal concave curvature radius. Default is
  \code{concave=convex}.
}
\item{resolution}{
  The internal computation resolution.  A warning will be issued when
  this needs to be increased for higher accuracy, with the required
  resolution stated.
}
\item{eps}{
  The polygonal curve simplification tolerance used for simplifying the
  resulting boundary curve.  See \code{\link{inla.simplify.curve}} for
  details.
}
}
\details{
  Morphological dilation by \code{convex}, followed by closing by
  \code{concave}, with minimum concave curvature radius \code{concave}.
  If the dilated set has no gaps of width between
   \deqn{2 convex (\sqrt{1+2 concave/convex} - 1)}{2*convex*(sqrt(1+2*concave/convex) - 1)}
 and
   \eqn{2 concave}{2*concave},
 then the minimum convex curvature radius is \code{convex}.
 Special case \code{concave=0} delegates to
 \code{inla.nonconvex.hull.basic}

 The implementation is based on the identity
   \deqn{dilation(a) \& closing(b) = dilation(a+b) \& erosion(b)}{dilation(a) & closing(b) = dilation(a+b) & erosion(b)}
 where all operations are with respect to disks with the specified radii.
  }
\value{
  An \code{\link{inla.mesh.segment}} object.
}
\author{
Finn Lindgren \email{finn.lindgren@gmail.com}
}
\note{
  Requires \code{nndistF} from the \code{splancs} package.
}

\examples{
if (require(splancs)) {
  loc = matrix(runif(20), 10, 2)
  boundary = inla.nonconvex.hull(loc, convex=0.2)
  lines(boundary, add=FALSE)
  points(loc)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
