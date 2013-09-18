\name{extract.groups}
\alias{extract.groups}
\alias{extract.groups.inla.mesh.segment}
\title{
Extract tagged boundary/internal segments.  
}
\description{
  Extract boundary or internal segments tagged by group id:s.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract.groups(...)

\method{extract.groups}{inla.mesh.segment}(
  segm, groups, groups.new = groups, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segm}{
    An \code{\link{inla.mesh.segment}} object.
%%     ~~Describe \code{\dots} here~~
}
  \item{groups}{
    The segment groups id:s to extract.
%%     ~~Describe \code{\dots} here~~
}
\item{groups.new}{
  Optional vector of group id remapping; \code{groups[k]} in the input
  will be replaced by \code{groups.new[k]} in the output.
%%     ~~Describe \code{\dots} here~~
}
\item{\dots}{
  Additional arguments, passed on to other methods.
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
\code{\link{inla.mesh.segment}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
