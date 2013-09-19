\name{inla.spde2.matern.sd.basis}
\alias{inla.spde2.matern.sd.basis}
\title{
  Approximate variance-compensating basis functions
}
\description{
  Calculates an approximate basis for \code{tau} and \code{kappa} for an
  \code{inla.spde2.matern} model where \code{tau} is a rescaling
  parameter.
}
\usage{
inla.spde2.matern.sd.basis(mesh, B.sd, B.range, method = 1,
                           local.offset.compensation = FALSE,
                           alpha = 2, ...)
}
\arguments{
  \item{mesh}{
    An \code{\link{inla.mesh}} object.
}
\item{B.sd}{
  Desired basis for log-standard deviations.
}
\item{B.range}{
  Desired basis for spatial range.
}
\item{method}{
  Construction method selector.  Expert option only.
}
\item{local.offset.compensation}{
  If \code{FALSE}, only compensate in the average for the \code{tau} offset.
}
\item{alpha}{
  The model \code{alpha} parameter.
}
\item{\dots}{
  Additional parameters passed on to internal \code{inla.spde2.matern}
  calls.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  List of basis specifications
  \item{B.tau }{Basis for \code{log(tau)}}
  \item{B.kappa }{Basis for \code{log(kappa)}}
  Intended for passing on to \code{\link{inla.spde2.matern}}.
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
\code{\link{inla.spde2.matern}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
