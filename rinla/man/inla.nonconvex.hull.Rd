\name{inla.nonconvex.hull}
\alias{inla.nonconvex.hull}
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
  The desired minimal concave curvature radius.
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
\value{
  An \code{\link{inla.mesh.segment}} object.
}
\author{
Finn Lindgren \email{finn.lindgren@gmail.com}
}
\note{
  Requires the \code{splancs} package.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
