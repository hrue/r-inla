\name{inla.spde2.generic}
\alias{inla.spde2.generic}
\alias{inla.spde2}
\alias{inla.spde2.theta2phi0}
\alias{inla.spde2.theta2phi1}
\alias{inla.spde2.theta2phi2}
\title{
  Generic spde2 model creation.
}
\description{
  Creates and inla.spde2 object describing the internal structure of an
  \code{'spde2'} model.
}
\usage{
inla.spde2.generic(M0, M1, M2, B0, B1, B2, theta.mu, theta.Q,
                   transform = c("logit", "log", "identity"),
                   theta.initial = theta.mu,
                   fixed = rep(FALSE, length(theta.mu)),
                   theta.fixed = theta.initial[fixed],
                   BLC = cbind(0, diag(nrow = length(theta.mu))),
                   ...)

## Map theta values to internal phi values
inla.spde2.theta2phi0(spde, theta)
inla.spde2.theta2phi1(spde, theta)
inla.spde2.theta2phi2(spde, theta)
}
\arguments{
  \item{M0}{
    The symmetric \code{M0} matrix.
}
  \item{M1}{
    The square \code{M1} matrix.
}
  \item{M2}{
    The symmetric \code{M2} matrix.
}
\item{B0}{
  Basis definition matrix for \eqn{\phi_0}{phi0}.
}
  \item{B1}{
  Basis definition matrix for \eqn{\phi_2}{phi1}.
}
  \item{B2}{
  Basis definition matrix for \eqn{\phi_2}{phi2}.
}
\item{theta.mu}{
  Prior expectation for the \eqn{\theta}{theta} vector
}
  \item{theta.Q}{
  Prior precision for the \eqn{\theta}{theta} vector
}
  \item{transform}{
    Transformation link for \eqn{\phi_2}{phi2}.  Valid settings are
    \code{"logit"}, \code{"log"}, and \code{"identity"}
}
  \item{theta.initial}{
  Initial value for the \eqn{\theta}{theta} vector. Default \code{theta.mu}
}
\item{fixed}{
  Logical vector.  For every \code{TRUE} value, treat the corresponding
  \code{theta} value as known.
}
\item{theta.fixed}{
  Vector holding the values of fixed \code{theta} values.  Default
  \code{=theta.initial[fixed]}
}
  \item{BLC}{
  Basis definition matrix for linear combinations of \code{theta}.
}
\item{\dots}{
  Additional parameters, currently unused.
}
\item{spde}{An \code{inla.sdpe2} object.}
\item{theta}{parameter values to be mapped.}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  For \code{inla.spde2.generic}, an \code{\link{inla.spde2}} object.

  For \code{inla.spde2.theta2phi0/1/2}, a vector of \eqn{\phi}{phi} values.
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
  \code{\link{inla.spde2.models}},
  \code{\link{inla.spde2.matern}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
