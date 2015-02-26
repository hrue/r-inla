\name{inla.fmesher.smorg}
\alias{inla.fmesher.smorg}
\title{
  Compute various mesh related quantities.
}
\description{
  Low level function for computing finite element matrices, spherical
  harmonics, B-splines, and point mappings with barycentric triangle
  coordinates.
}
\usage{
inla.fmesher.smorg(loc, tv,
                   fem = NULL, aniso = NULL,
                   gradients = FALSE,
                   sph0 = NULL, sph = NULL, bspline = NULL,
                   points2mesh = NULL,
                   splitlines = NULL,
                   output = NULL,
                   keep = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc}{
    3-column triangle vertex coordinate matrix.
}
\item{tv}{
  3-column triangle vertex index matrix.
}
\item{fem}{
  Maximum finite element matrix order to be computed.
}
\item{aniso}{
  A two-element list with \eqn{\gamma}{gamma} and \eqn{v}{v} for an
  anisotropic operator \eqn{\nabla\cdot H \nabla}{div H grad}, where
  \eqn{H=\gamma I + v v^\top}{H = gamma I + v v'}
}
\item{gradients}{
  When \code{TRUE}, calculate derivative operator matrices \code{dx},
  {dy}, and {dz}.
}
\item{sph0}{
  Maximal order of rotationally invariant spherical harmonics.
}
  \item{sph}{
  Maximal order of general spherical harmonics.
}
\item{bspline}{
  Rotationally invariant B-splines on a sphere.
  3-vector with number of basis functions \code{n}, basis degree
  \code{degree}, and a logical; \code{TRUE} uniform knot angles,
  \code{FALSE} for uniform spacing in \eqn{\sin(latitude)}{sin(latitude)}.
}
\item{points2mesh}{
  3-column matrix with points to be located in the mesh.
}
\item{splitlines}{
  A list with elements \code{loc} (3-column coordinate matrix) and
  \code{idx} (2-column index matrix) describing line segments that are
  to be split into sub-segments at triangle boundaries.
}
\item{output}{
  Names of objects to be included in the output, if different from defaults.
}
\item{keep}{
  When \code{TRUE}, for debugging purposes keep the \code{fmesher} I/O
  files on disk.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of generated named quantities.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ internal }
