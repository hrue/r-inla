\name{inla.mesh.deriv}
\alias{inla.mesh.deriv}
\title{
  Directional derivative matrices for functions on meshes.
}
\description{
  Calculates directional derivative matrices for functions on
  \code{\link{inla.mesh}} objects.
}
\usage{
inla.mesh.deriv(mesh, loc)
}
\arguments{
  \item{mesh}{
    An \code{\link{inla.mesh}} object.
}
  \item{loc}{
    Coordinates where the derivatives should be evaluated.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{A }{The projection matrix, \code{u(loc_i)=sum_j A_ij w_i}}
\item{dx, dy, dz }{Derivative weight matrices, \code{du/dx(loc_i)=sum_j
    dx_ij w_i}, etc.}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
