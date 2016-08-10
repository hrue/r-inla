\name{inla.CRS}
\alias{inla.CRS}
\title{
  Create a coordinate reference system object
}
\description{
  Creates either a CRS object or an inla.CRS object, describing a
  coordinate reference system 
}
\usage{
inla.CRS(projargs = NA_character_, doCheckCRSArgs = TRUE,
         args = NULL, oblique = NULL, ...)
}
\arguments{
  \item{projargs}{
    Either 1) a projection argument string suitable as input to
    \code{sp::CRS}, or 2) an existing \code{CRS} object,
    or 3) a shortcut reference string to a predefined projection
    (\code{longlat}, \code{lambert}, \code{mollweide},
    \code{hammer}, and \code{sphere}).
}
\item{doCheckCRSArgs}{
  default TRUE, must be set to FALSE by package
  developers including \code{CRS} in an S4 class definition to avoid
  uncontrolable loading of the \code{rgdal} namespace.
}
\item{args}{
  An optional list of name/value pairs to add to and/or override the
  PROJ4 arguments in \code{projargs}.  \code{name=value} is converted to
  \code{"+name=value"}, and \code{name=NA} is converted to
  \code{"+name"}.
}
\item{oblique}{
  Vector of length at most 4 of rotation angles (in degrees) for an oblique
  projection, all values defaulting to zero. The values indicate
  (longitude, latitude, orientation, orbit), as explained in the Details
  section below.
}
\item{\dots}{Additional parameters. Not currently in use.}
}
\details{
  %  ~~ If necessary, more details than the description above ~~
  The first two elements of the \code{oblique} vector are
  the (longitude, latitude) coordinates for the oblique centre
  point. The third value (orientation) is a counterclockwise rotation
  angle for an observer looking at the centre point from outside the
  sphere. The fourth value is the quasi-longitude (orbit angle) for
  a rotation along the oblique observers equator.

  Simple oblique: \code{oblique=c(0, 45)}

  Polar: \code{oblique=c(0, 90)}

  Quasi-transversal: \code{oblique=c(0, 0, 90)}

  Satellite orbit viewpoint:
  \code{oblique=c(lon0-time*v1, 0, orbitangle, orbit0+time*v2)},
  where \code{lon0} is the longitude at which a
  satellite orbit crosses the equator at \code{time=0}, when the
  satellite is at an angle \code{orbit0} further along in its orbit.
  The orbital angle relative to the equatorial plane is
  \code{orbitangle}, and \code{v1} and \code{v2} are the angular
  velocities of the planet and the satellite, respectively.
  Note that "forward" from the satellite's point of view is
  "to the right" in the projection.

  When \code{oblique[2]} or
  \code{oblique[3]} are non-zero, the resulting projection is only
  correct for perfect spheres.
}
\value{
  Either an \code{sp::CRS} object or an \code{inla.CRS} object,
  depending on if the coordinate reference system described by the
  parameters can be expressed with a pure \code{sp::CRS} object or not.

  An S3 \code{inla.CRS} object is a list, usually (but not necessarily)
  containing at least one element:
  \item{crs }{The basic \code{sp::CRS} object}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Finn Lindgren <finn.lindgren@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[sp]{CRS}}, \code{\link{inla.CRSargs}},
  \code{\link{plot.CRS}}, \code{\link{inla.identical.CRS}}
}
\examples{
if (require(rgdal)) {
  halfroot <- "+a=0.7071067811865476 +b=0.7071067811865476"
  crs1 <- inla.CRS("+proj=longlat +ellps=sphere +a=1 +b=1")
  crs2 <- inla.CRS("+proj=cea +ellps=sphere +lat_ts=0 +units=m +a=1 +b=1")
  crs3 <- inla.CRS(paste("+proj=moll +ellps=sphere +units=m", halfroot))
  crs4 <- inla.CRS(paste("+proj=hammer +ellps=sphere +units=m", halfroot))
  crs5 <- inla.CRS("+proj=geocent +ellps=sphere +a=1 +b=1 +units=m")
  ## Shortcuts:
  crs1 <- inla.CRS("longlat")
  crs2 <- inla.CRS("lambert")
  crs3 <- inla.CRS("mollweide")
  crs4 <- inla.CRS("hammer")
  crs5 <- inla.CRS("sphere")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
