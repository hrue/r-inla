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
         params = NULL, orient = NULL)
}
\arguments{
  \item{projargs}{
    Either a projection argument string suitable as input to
    \code{sp::CRS}, or a shortcut reference string to a predefined projection
    (\code{longlat}, \code{sphere}, \code{mollweide}, \code{lambert}).
}
\item{doCheckCRSArgs}{
  default TRUE, must be set to FALSE by package
  developers including \code{CRS} in an S4 class definition to avoid
  uncontrolable loading of the \code{rgdal} namespace.
}
\item{params}{
  Projection parameters. Not currently in use.
}
\item{orient}{
  Vector of length 4 of rotation angles for an oblique projection.
  Not currently in use.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Either an \code{sp::CRS} object or an \code{inla.CRS} object,
  depending on if the coordinate reference system described by the
  parameters can be expressed with a pure \code{sp::CRS} object or not.

  An S3 \code{inla.CRS} object is a list containing at least one element: 
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
  \code{\link{sp:CRS}}, \code{\link{inla.CRSargs}}
}
\examples{
crs1 <- inla.CRS("+proj=longlat +ellps=sphere +a=1 +b=1")
crs2 <- inla.CRS("+proj=geocent +ellps=sphere +a=1 +b=1 +units=m", 
crs3 <- inla.CRS("+proj=moll +ellps=sphere +units=m +a=0.7071067811865476 +b=0.7071067811865476", 
crs4 <- inla.CRS("+proj=cea +ellps=sphere +lat_ts=0 +units=m +a=1 +b=1")
## Shortcuts:
crs1 <- inla.CRS("longlat"), 
crs2 <- inla.CRS("sphere"), 
crs3 <- inla.CRS("mollwide"), 
crs4 <- inla.CRS("lambert"), 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
