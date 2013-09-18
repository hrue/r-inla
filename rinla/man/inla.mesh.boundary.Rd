\name{inla.mesh.boundary}
\alias{inla.mesh.boundary}
\alias{inla.mesh.interior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Constraint segment extraction for inla.mesh
}
\description{
  Constructs an list of \code{inla.mesh.segment} object from boundary or
  interior constraint information in an \code{\link{inla.mesh}} object.
}
\usage{
inla.mesh.boundary(mesh, grp = NULL)
inla.mesh.interior(mesh, grp = NULL)
}
\arguments{
  \item{mesh}{
    An \code{inla.mesh} object.
  }
  \item{grp}{
    Group indices to extract.  If \code{NULL}, all boundary/interior
    constrain groups are extracted.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of \code{inla.mesh.segment} objects.
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
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh.segment}},
  \code{\link{inla.mesh.create}},
  \code{\link{inla.mesh.create.helper}}
}
\examples{
loc = matrix(runif(100*2)*1000,100,2)
mesh = inla.mesh.create.helper(points.domain=loc, max.edge=c(50,500))
boundary = inla.mesh.boundary(mesh)
interior = inla.mesh.interior(mesh)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
