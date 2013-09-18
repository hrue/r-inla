\name{inla.mesh.map}
\alias{inla.mesh.map}
\alias{inla.mesh.map.lim}
\title{
  Coordinate mappings for \code{inla.mesh} projections.
}
\description{
  Calculates coordinate mappings for \code{inla.mesh} projections.
}
\usage{
inla.mesh.map(loc,
              projection = c("default", "longlat",
                             "longsinlat", "mollweide"),
              inverse = TRUE)

## Compute sensible default map axis limits
inla.mesh.map.lim(loc = NULL,
                  projection = c("default", "longlat",
                                 "longsinlat", "mollweide"))
}
\arguments{
  \item{loc}{
    Coordinates to be mapped.
}
\item{projection}{
  The projection type.
}
\item{inverse}{
  If \code{TRUE}, \code{loc} are map coordinates and coordinates in the
  mesh domain are calculated.
  If \code{FALSE}, \code{loc} are coordinates in the mesh domain and the
  forward map projection is calculated.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  For \code{inla.mesh.map.lim}, a list:
  \item{xlim }{X axis limits in the map domain}
  \item{ylim }{Y axis limits in the map domain}
  No attempt is made to find minimal limits for partial spherical domains.
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

\seealso{
  \code{\link{inla.mesh.project}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
