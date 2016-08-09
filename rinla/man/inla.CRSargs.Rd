\name{inla.CRSargs}
\alias{inla.CRSargs}
\alias{inla.as.list.CRS}
\alias{inla.as.list.CRSargs}
\alias{inla.as.CRS.list}
\alias{inla.as.CRSargs.list}
\title{
  Show expanded CRS arguments
}
\description{
Wrapper for \code{sp::CRS} and \code{inla.CRS} objects to extract the 
coordinate reference system argument string.  
}
\usage{
inla.CRSargs(x, ...)
inla.as.list.CRS(x, ...)
inla.as.list.CRSargs(x, ...)
inla.as.CRS.list(x, ...)
inla.as.CRSargs.list(x, ...)
}
\arguments{
  \item{x}{
    An \code{sp::CRS} or \code{inla.CRS} object (for
    \code{inla.CRSargs} and \code{inla.as.list.CRS}), a character string
    (for \code{inla.as.list.CRSargs}), or a list (for
    \code{inla.as.CRS.list} and \code{inla.as.CRSargs.list}).
  }
  \item{\dots}{Additional arguments passed on to other methods.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  For \code{inla.CRSargs} and \code{inla.as.CRSargs.list}, a character string with PROJ.4 arguments.

  For \code{inla.as.list.CRS} and \code{inla.as.list.CRSargs}, a list of
  name/value pairs.

  For \code{inla.as.CRS.list}, a \code{CRS} or \code{inla.CRS} object.
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
\code{\link[rgdal]{CRSargs}},
\code{\link{inla.CRS}}
}
\examples{
if (require(rgdal)) {
  crs0 <- inla.CRS("longlat")
  p4s <- inla.CRSargs(crs0)
  lst <- inla.as.list.CRSargs(p4s)
  crs1 <- inla.as.CRS.list(lst)
  lst$a <- 2
  crs2 <- inla.CRS(p4s, args=lst)
  print(inla.CRSargs(crs0))
  print(inla.CRSargs(crs1))
  print(inla.CRSargs(crs2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
