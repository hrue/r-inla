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
     outline = TRUE,
     graticule = c(15, 15, 45),
     tissot = c(30, 30, 30),
     asp = 1, add = FALSE,
     eps=0.05, ...)
\method{plot}{CRS}(x,
     xlim = NULL, ylim = NULL,
     outline = TRUE,
     graticule = c(15, 15, 45),
     tissot = c(30, 30, 30),
     asp = 1, add = FALSE,
     eps=0.05, ...)
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
  Vector of length at most 3, to plot meridians with spacing
  \code{graticule[1]} degrees and parallels with spacing
  \code{graticule[2]} degrees. \code{graticule[3]} optionally specifies
  the spacing above and below the first and last parallel.
  When \code{graticule[1]==0} no meridians
  are drawn, and when \code{graticule[2]==0} no parallels are drawn. Use
  \code{graticule=NULL} to skip drawing a graticule.
}
\item{tissot}{
  Vector of length at most 3, to plot Tissot's indicatrices with spacing
  \code{tissot[1]} degrees and parallels with spacing
  \code{tissot[2]} degrees. \code{tissot[3]} specifices a scaling
  factor.
  Use \code{tissot=NULL} to skip drawing a Tissot's indicatrices.
}
\item{asp}{
  The aspect ratio for the plot, default 1.
}
\item{add}{
  If \code{TRUE}, add the projecton plot to an existing plot.
}
\item{eps}{
  Clipping tolerance for rudimentary boundary clipping
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
  oblique <- c(0,45,45,0)
  for (projtype in c("longlat", "lambert", "mollweide", "hammer")) {
    plot(inla.CRS(projtype), main=projtype)
    plot(inla.CRS(projtype, oblique=oblique), main=paste("oblique", projtype))
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
