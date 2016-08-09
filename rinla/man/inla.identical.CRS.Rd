\name{inla.identical.CRS}
\alias{inla.identical.CRS}
\title{
  Test CRS and inla.CRS for equality
}
\description{
  Wrapper for identical, optionally testing only the CRS part of two objects
}
\usage{
inla.identical.CRS(crs0, crs1, crsonly = FALSE)
}
\arguments{
  \item{crs0}{
    A \code{CRS} or \code{inla.CRS} object.
}
  \item{crs1}{
    A \code{CRS} or \code{inla.CRS} object.
}
\item{crsonly}{
  Logical. If \code{TRUE}, only the \code{CRS} part of a \code{inla.CRS}
  object is compared.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
crs0 <- inla.CRS("longlat")
crs1 <- inla.CRS("longlat", oblique=c(0,90))
print(c(inla.identical.CRS(crs0, crs0),
        inla.identical.CRS(crs0, crs1),
        inla.identical.CRS(crs0, crs1, crsonly=TRUE)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
