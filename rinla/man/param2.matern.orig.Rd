\name{param2.matern.orig}
\alias{param2.matern.orig}
\title{
  Parameter settings for \code{inla.spde2.matern} models.
}
\description{
  Construct parameter settings for \code{inla.spde2.matern} models.
}
\usage{
param2.matern.orig(mesh, alpha = 2,
                   B.tau = matrix(c(0, 1, 0), 1, 3),
                   B.kappa = matrix(c(0, 0, 1), 1, 3),
                   prior.variance.nominal = 1,
                   prior.range.nominal = NULL,
                   prior.tau = NULL,
                   prior.kappa = NULL,
                   theta.prior.mean = NULL,
                   theta.prior.prec = 0.1)
}
\arguments{
  \item{mesh}{
    The mesh to build the model on, as an \code{\link{inla.mesh}} object.
  }
  \item{alpha}{Fractional operator order, \eqn{0<\alpha\leq 2}
    supported. (\eqn{\nu=\alpha-d/2})}
  \item{B.tau}{Matrix with specification of log-linear model for
    \eqn{\tau}.}
  \item{B.kappa}{Matrix with specification of log-linear model for
    \eqn{\kappa}.}
  \item{prior.variance.nominal}{Nominal prior mean for the field variance}
  \item{prior.range.nominal}{Nominal prior mean for the spatial range}
  \item{prior.tau}{Prior mean for tau (overrides \code{prior.variance.nominal})}
  \item{prior.kappa}{Prior mean for kappa (overrides \code{prior.range.nominal})}
  \item{theta.prior.mean}{(overrides \code{prior.*})}
  \item{theta.prior.prec}{Scalar, vector or matrix, specifying the joint
    prior precision for \eqn{theta}.}
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
\code{\link{inla.spde2.matern}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
