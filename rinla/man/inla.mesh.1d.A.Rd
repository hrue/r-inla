\name{inla.mesh.1d.A}
\alias{inla.mesh.1d.A}
\alias{inla.mesh.1d.bary}
\title{
  Mapping matrix for 1D meshes
}
\description{
  Calculates barycentric coordinates and weight matrices for
  \code{\link{inla.mesh.1d}} objects.
}
\usage{
inla.mesh.1d.A(mesh, loc, weights = NULL, derivatives = NULL,
               method = c("linear", "nearest", "quadratic"))

inla.mesh.1d.bary(mesh, loc, method = c("linear", "nearest"))
}
\arguments{
  \item{mesh}{
    An \code{\link{inla.mesh.1d}} object.
}
\item{loc}{
  Coordinate values.
}
\item{weights}{
  Weights to be applied to the \code{A} matrix rows.
}
\item{derivatives}{
  If \code{TRUE}, also compute derivative weight matrices \code{dA} and
  \code{d2A}.
}
\item{method}{
  Interpolation method.  If not specified for \code{inla.mesh.1d.A}
  (recommended), it is determined by the \code{mesh} basis function
  properties.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
