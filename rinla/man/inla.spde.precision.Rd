\name{inla.spde.precision}
\alias{inla.spde.precision}
\alias{inla.spde1.precision}
\alias{inla.spde2.precision}
\alias{inla.spde.precision.inla.spde1}
\alias{inla.spde.precision.inla.spde2}
\title{
  Precision matrices for SPDE models
}
\description{
Calculates the precision matrix for given parameter values based on an
\code{inla.spde} model object.
}
\usage{
%% Standard usage for inla.spde2.models:
%% inla.spde.precision(spde, theta = c(...))

inla.spde.precision(...)

\method{inla.spde.precision}{inla.spde2}(spde,
                    theta = NULL,
                    phi0 = inla.spde2.theta2phi0(spde, theta),
                    phi1 = inla.spde2.theta2phi1(spde, theta),
                    phi2 = inla.spde2.theta2phi2(spde, theta), ...)
inla.spde2.precision(spde,
                     theta = NULL,
                     phi0 = inla.spde2.theta2phi0(spde, theta),
                     phi1 = inla.spde2.theta2phi1(spde, theta),
                     phi2 = inla.spde2.theta2phi2(spde, theta), ...)

## For deprecated inla.spde1 models:
\method{inla.spde.precision}{inla.spde1}(spde, ...)
inla.spde1.precision(spde, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spde}{
    An \code{inla.spde} object.
  }
  \item{theta}{
    The parameter vector.
  }
  \item{phi0}{
    Internal parameter for a generic model.  Expert option only.
  }
  \item{phi1}{
    Internal parameter for a generic model.  Expert option only.
  }
  \item{phi2}{
    Internal parameter for a generic model.  Expert option only.
  }
  \item{\dots}{
    Additional parameters passed on to other methods.
    }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A sparse precision matrix.
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
\code{\link{inla.spde.models}},
\code{\link{inla.spde2.generic}},
\code{\link{inla.spde2.theta2phi0}},
\code{\link{inla.spde2.theta2phi1}},
\code{\link{inla.spde2.theta2phi2}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
