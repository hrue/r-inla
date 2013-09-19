\name{inla.spde.sample}
\alias{inla.spde.sample}
\alias{inla.spde.sample.default}
\alias{inla.spde.sample.inla.spde}
\title{
  Sample from SPDE models
}
\description{
  Old methods fo sampling from a SPDE model.  For new code, use
  \code{\link{inla.spde.precision}} and \code{\link{inla.qsample}}
  instead.
}
\usage{
inla.spde.sample(...)

\method{inla.spde.sample}{default}(precision, seed=NULL, ...)
\method{inla.spde.sample}{inla.spde}(spde, seed=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{precision}{
    A precision matrix.
  }
  \item{seed}{
    The seed for the pseudo-random generator.
  }
  \item{spde}{
    An \code{inla.spde} object.
  }
  \item{\dots}{
    Parameters passed on to other methods.
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

\seealso{
  \code{\link{inla.spde.precision}},
  \code{\link{inla.qsample}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
