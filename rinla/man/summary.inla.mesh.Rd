\name{summary.inla.mesh}
\alias{summary.inla.mesh}
\alias{print.summary.inla.mesh}
\title{
  Summarizing triangular mesh objects
}
\description{
  Construct and print \code{inla.mesh} object summaries
}
\usage{
\method{summary}{inla.mesh}(object, verbose = FALSE, ...)

\method{print}{summary.inla.mesh}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class \code{"inla.mesh"}, usually a result of a call to
  \code{\link{inla.mesh.create}} or \code{\link{inla.mesh.2d}}.
}
  \item{x}{
  an object of class \code{"summary.inla.mesh"}, usually a result of a call to
  \code{\link{summary.inla.mesh}}.
}
\item{verbose}{
  If \code{TRUE}, produce a more detailed output. 
}
\item{\dots}{
  further arguments passed to or from other methods.
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
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
