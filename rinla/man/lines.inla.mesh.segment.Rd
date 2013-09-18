\name{lines.inla.mesh.segment}
\alias{lines.inla.mesh.segment}
\title{
  Draw \code{inla.mesh.segment} objects.
}
\description{
  Draws a \code{\link{inla.mesh.segment}} object with generic or
  \code{rgl} graphics.
}
\usage{
\method{lines}{inla.mesh.segment}(x, loc = NULL, col = NULL,
     colors = c("black", "blue", "red", "green"),
     add = TRUE, xlim = NULL, ylim = NULL,
     rgl = FALSE, ...)
}
\arguments{
  \item{x}{
    An \code{\link{inla.mesh.segment}} object.
}
\item{loc}{
  Point locations to be used if \code{x$loc} is \code{NULL}.
}
\item{col}{
  Segment color specification.
}
\item{colors}{
  Colors to cycle through if \code{col} is \code{NULL}.
}
\item{add}{
  If \code{TRUE}, add to the current plot, otherwise start a new plot.
}
\item{xlim}{
  X axis limits for a new plot.
}
  \item{ylim}{
  Y axis limits for a new plot.
}
\item{rgl}{
  If \code{TRUE}, use \code{rgl} for plotting.
}
\item{\dots}{
  Additional parameters, passed on to graphics methods.
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
Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh.segment}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
