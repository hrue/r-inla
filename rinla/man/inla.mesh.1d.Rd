\name{inla.mesh.1d}
\alias{inla.mesh.1d}
\title{
  Function space definition objects for 1D SPDE models.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inla.mesh.1d(loc,
             interval = range(loc),
             boundary = NULL,
             degree = 1,
             free.clamped = FALSE,
             cyclic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc}{
    B-spline knot locations.
}
\item{interval}{
  Interval domain endpoints.
}
\item{boundary}{
  Boundary condition specification.  Valid conditions are
  \code{c('neumann', 'dirichlet', 'free', 'cyclic')}.  Two separate
  values can be specified, one applied to each endpoint.
}
\item{degree}{
  The B-spline basis degree.  Supported values are 0, 1, and 2.
}
\item{free.clamped}{
  If \code{TRUE}, for \code{'free'} boundaries, clamp the basis
  functions to the interval endpoints.
}
\item{cyclic}{
  Deprecated option, kept for backwards compatibility.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Finn Lindgren \email{finn.lindgren@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
