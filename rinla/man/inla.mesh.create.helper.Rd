\name{inla.mesh.create.helper}
\alias{inla.mesh.create.helper}
\title{
High-quality triangulations
}
\description{
  Create a triangle mesh based on initial point locations, specified or
  automatic boundaries, and mesh quality parameters.
}
\usage{
inla.mesh.create.helper(points = NULL,
                        points.domain = NULL,
                        ...)
}
\arguments{
  \item{points}{
    Matrix of point locations to be used as initial triangulation nodes.
  }
  \item{points.domain}{
    Matrix of point locations used to determine the domain extent.
  }
  \item{\dots}{
    Additional arguments passed on to \code{\link{inla.mesh.2d}}.
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
\note{
  Obsolete function kept for backwards compatibility.  Use
  \code{\link{inla.mesh.2d}} instead.
}

\seealso{
\code{\link{inla.mesh.2d}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
