\name{inla.simplify.curve}
\alias{inla.simplify.curve}
\title{
  Recursive curve simplification.
}
\description{
  Attempts to simplify a polygonal curve by joining nearly colinear segments.
}
\usage{
inla.simplify.curve(loc, idx, eps)
}
\arguments{
  \item{loc}{
    Coordinate matrix.
}
\item{idx}{
  Index vector into \code{loc} specifying a polygonal curve.
}
\item{eps}{
  Straightness tolerance.
}
}
\details{
  Uses a variation of the binary splitting Ramer-Douglas-Peucker
  algorithm, with a width \code{eps} ellipse instead of a rectangle,
  motivated by prediction ellipse for Brownian bridge.
}
\value{
  An index vector into \code{loc} specifying the simplified polygonal curve.
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
\examples{
theta = seq(0, 2*pi, length=1000)
loc = cbind(cos(theta), sin(theta))
idx = inla.simplify.curve(loc=loc, idx=1:nrow(loc), eps=0.01)
print(c(nrow(loc), length(idx)))
plot(loc, type="l")
lines(loc[idx,], col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
