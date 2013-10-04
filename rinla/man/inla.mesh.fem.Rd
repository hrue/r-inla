\name{inla.mesh.fem}
\alias{inla.mesh.fem}
\alias{inla.mesh.1d.fem}
\title{
  Finite element matrices
}
\description{
Constructs finite element matrices for \code{\link{inla.mesh}} and
\code{\link{inla.mesh.1d}} objects.
}
\usage{
## 2D and 1D meshes
inla.mesh.fem(mesh, order = 2)

## 1D meshes, order 2 models only
inla.mesh.1d.fem(mesh)
}
\arguments{
  \item{mesh}{
    An \code{\link{inla.mesh}} or \code{\link{inla.mesh.1d}} object.
}
\item{order}{
  The model order.
}
}
%\details{  
%}
\value{
  A list of sparse matrices based on basis functions \code{psi_i}:
  \item{c0 }{\code{c0[i,j] = < psi_i, 1 >}}
  \item{c1 }{\code{c1[i,j] = < psi_i, psi_j >}}
  \item{g1 }{\code{g1[i,j] = < grad psi_i, grad psi_j >}}
  \item{g2 }{\code{g2 = g1 * c0^-1 * g1}}
  \item{gk }{\code{gk = g1 * (c0^-1 * g1)^(k-1)}, up to and including
    \code{k=order}}
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
