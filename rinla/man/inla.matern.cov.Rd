\name{inla.matern.cov}
\alias{inla.matern.cov}
\alias{inla.matern.cov.s2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Numerical evaluation of Matern and related covariance functions.
}
\description{
  Calculates covariance and correlation functions for Matern models
  and related oscillating SPDE models, on \eqn{R^d}{R^d} and on the
  sphere, \eqn{S^2}{S^2}.
}
\usage{
inla.matern.cov(nu, kappa, x,
                d = 1,
                corr = FALSE,
                norm.corr = FALSE,
                theta,
                epsilon = 1e-08)

inla.matern.cov.s2(nu, kappa, x,
                   norm.corr = FALSE,
                   theta = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{
    The Matern smoothness parameter.
  }
  \item{kappa}{
    The spatial scale parameter.
  }
  \item{x}{
    Distance values.
  }
  \item{d}{
    Space dimension; the domain is \eqn{R^d}{R^d}.
  }
  \item{corr}{
    If \code{TRUE}, calculate correlations, otherwise calculate
    covariances.  Only used for pure Matern models
    (i.e. with \eqn{\theta=0}{theta=0}).
  }
  \item{norm.corr}{
    If \code{TRUE}, normalise by the estimated variance, giving
    approximate correlations.
  }
  \item{theta}{
    Oscillation strength parameter.
  }
  \item{epsilon}{
    Tolerance for detecting points close to distance zero.
  }
}
\details{
  On \eqn{R^d}{R^d}, the models are \emph{defined} by the spectral density
  given by
  \deqn{S(w) = \frac{1}{(2\pi)^d (\kappa^4 + 2 \kappa^2 \cos(\pi  
      \theta) |w|^2 + |w|^4)^{(\nu + d/2)/2}}
    }{S(w) = 1 / ( (2\pi)^d * (kappa^4 + 2 kappa^2 * cos(pi * theta) *
      |w|^2 + |w|^4)^((nu + d/2)/2) )}
    
    On \eqn{S^2}{S^2}, the models are \emph{defined} by the spectral
  coefficients 
  \deqn{S(k) = \frac{2k+1}{4\pi (\kappa^4 + 2 \kappa^2 \cos(\pi  
      \theta) k(k+1) + k^2(k+1)^2)^{(\nu + 1)/2}}
      }{S(k) = (2k+1) / (4 pi (kappa^4 + 2 kappa^2 cos(pi  
      theta) k(k+1) + k^2(k+1)^2)^((\nu + 1)/2) )}
}
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
%%\keyword{ ~kwd1 }
