\name{inla.mesh}
\alias{inla.mesh}
\title{
  Create high-quality triangulations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inla.mesh(loc = NULL,
          boundary = NULL, interior = NULL,
          extend = (missing(tv) || is.null(tv)),
          refine = FALSE,
          tv = NULL, grid = NULL, manifold = NULL,
          cutoff.distance = 0,
          plot.delay = NULL,
          dir = tempdir(), prefix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc}{
}
  \item{boundary}{
%%     ~~Describe \code{boundary} here~~
}
  \item{interior}{
%%     ~~Describe \code{interior} here~~
}
  \item{extend}{
%%     ~~Describe \code{extend} here~~
}
  \item{refine}{
%%     ~~Describe \code{refine} here~~
}
  \item{tv}{
%%     ~~Describe \code{tv} here~~
}
  \item{grid}{
%%     ~~Describe \code{grid} here~~
}
  \item{manifold}{
%%     ~~Describe \code{manifold} here~~
}
  \item{cutoff.distance}{
%%     ~~Describe \code{cutoff.distance} here~~
}
  \item{plot.delay}{
%%     ~~Describe \code{plot.delay} here~~
}
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An \code{fmesher.mesh} object.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{mesh.segm}}, \code{\link{mesh.grid}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fmesher }
