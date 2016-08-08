\name{inla.spTransform}
\alias{inla.spTransform}
\alias{inla.spTransform.default}
\alias{inla.spTransform.SpatialPoints}
\alias{inla.spTransform.inla.mesh.lattice}
\alias{inla.spTransform.inla.mesh.segment}
\alias{inla.spTransform.inla.mesh}
\title{
  Wrapper method for \code{sp::spTransform}
}
\description{
  Handles transformation of various inla objects accorting to coordinate
  reference systems of \code{sp::CRS} or \code{inla.CRS} class.
}
\usage{
inla.spTransform(x, ...)
\method{inla.spTransform}{default}(x, crs0, crs1,
                 passthrough=FALSE, ...)
\method{inla.spTransform}{SpatialPoints}(x, CRSobj,
                 passthrough=FALSE, ...)
\method{inla.spTransform}{inla.mesh.lattice}(x, CRSobj,
                 passthrough=FALSE, ...)
\method{inla.spTransform}{inla.mesh.segment}(x, CRSobj,
                 passthrough=FALSE, ...)
\method{inla.spTransform}{inla.mesh}(x, CRSobj,
                 passthrough=FALSE, ...)
}
\arguments{
  \item{x}{
    The object that should be transformed from it's current CRS to a new CRS
}
\item{crs0}{
  The source \code{sp::CRS} or \code{inla.CRS} object
}
  \item{crs1}{
  The target \code{sp::CRS} or \code{inla.CRS} object
}
  \item{CRSobj}{
  The target \code{sp::CRS} or \code{inla.CRS} object
}
\item{passthrough}{
  default FALSE. Setting to TRUE allows objects with no CRS information
  to be passed through without transformation.
}
\item{\dots}{
  Potential additional arguments
}
}
%\details{ 
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  The object is returned with its coordinates transformed
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
  Finn Lindgren <finn.lindgren@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.CRS}}
}
\examples{
if (require(rgdal)) {
  latt <- inla.mesh.lattice(-10:10, 40:60)
  mesh1 <- inla.mesh.create(lattice=latt, extend=FALSE, refine=FALSE,
                            crs=inla.CRS("longlat"))
  mesh2 <- inla.spTransform(mesh1, inla.CRS("lambert"))
  summary(mesh1)
  summary(mesh2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
