\name{plot.inla.CRS}
\alias{plot.inla.CRS}
\alias{plot.CRS}
\title{
  Plot CRS and inla.CRS objects
}
\description{
  Plot the outline of a CRS or inla.CRS projection, with optional
  graticules (transformed parallels and meridians) and Tissot indicatrices.
}
\usage{
\method{plot}{inla.CRS}(x,
     xlim = NULL, ylim = NULL,
     outline = TRUE, graticule = c(24, 12), tissot = c(12, 6),
     asp = 1, add = FALSE, ...)
\method{plot}{CRS}(x, ...)
}
\arguments{
  \item{x}{
    A \code{CRS} or \code{\link{inla.CRS}} object.
}
\item{xlim}{
  Optional x-axis limits.
}
  \item{ylim}{
  Optional y-axis limits.
}
\item{outline}{
  Logical, if \code{TRUE}, draw the outline of the projection.
}
\item{graticule}{
  If a length 2 vector is given, plot meridians (spacing
  \code{360/n[1]}) and parallels (spacing \code{180/n[2]}). Set
  \code{graticule=NULL} to skip drawing graticules.
}
\item{tissot}{
  Not yet implemented.
}
\item{asp}{
  The aspect ratio for the plot, default 1.
}
\item{add}{
  If \code{TRUE}, add the projecton plot to an existing plot.
}
\item{\dots}{
  Additional arguments passed on to the internal calls to \code{plot}
  and \code{lines}.
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
if (require(rgdal)) {
  orient <- c(0,45,45,0)
  par(mfrow=c(4,2))
  for (projtype in c("longlat", "lambert", "mollweide", "hammer")) {
    plot(inla.CRS(projtype), main=projtype)
    plot(inla.CRS(projtype, orient=orient), main=paste("oblique", projtype))
  }
  par(mfrow=c(1,1))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
