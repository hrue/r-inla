\name{inla.diameter}
\alias{inla.diameter}
\alias{inla.diameter.default}
\alias{inla.diameter.inla.mesh}
\alias{inla.diameter.inla.mesh.segment}
\alias{inla.diameter.inla.mesh.lattice}
\alias{inla.diameter.inla.mesh.1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Diameter of a point set
}
\description{
  Find an upper bound to the convex hull of a point set
}
\usage{
inla.diameter(x, ...)

\method{inla.diameter}{default}(x, manifold="", ...)
\method{inla.diameter}{inla.mesh}(x, ...)
\method{inla.diameter}{inla.mesh.segment}(x, ...)
\method{inla.diameter}{inla.mesh.lattice}(x, ...)
\method{inla.diameter}{inla.mesh.1d}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A point set as an \eqn{n\times d}{n x d} matrix, or an
    \code{inla.mesh} related object.
}
\item{manifold}{
  Character string specifying the manifold type. Default is to treat the
  point set with Euclidean \eqn{R^d} metrics. Use \code{manifold="S2"}
  for great circle distances on the unit sphere (this is set
  automatically for \code{inla.mesh} objects).
}
\item{\dots}{
  Additional parameters passed on to other methods.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A scalar, upper bound for the diameter of the convex hull of the point set.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Finn Lindgren <finn.lindgren@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
inla.diameter(matrix(c(0,1,1,0, 0,0,1,1), 4, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
