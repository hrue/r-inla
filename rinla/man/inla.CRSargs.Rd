\name{inla.CRSargs}
\alias{inla.CRSargs}
\title{
  Show expanded CRS arguments
}
\description{
Wrapper for \code{sp::CRS} and \code{inla.CRS} objects to extract the 
coordinate reference system argument string.  
}
\usage{
inla.CRSargs(crs)
}
\arguments{
  \item{crs}{
    An \code{sp::CRS} or \code{inla.CRS} object
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A character string with PROJ.4 arguments.
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
\code{\link{rgdal::CRSargs}},
\code{\link{inla.CRS}}
}
\examples{
inla.CRSargs(inla.CRS("longlat"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
