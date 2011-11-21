\name{inla.spde.make.A}
\alias{inla.spde.make.A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add title
%%  ~~function to do ... ~~
}
\description{Add description
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inla.spde.make.A(mesh = NULL,
                 loc = NULL,
                 index = NULL,
                 group = NULL,
                 repl = 1L,
                 n.mesh = NULL,
                 n.group = max(group),
                 n.repl = max(repl),
                 group.mesh = NULL,
                 group.method = c("nearest", "S0", "S1"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{
%%     ~~Describe \code{mesh} here~~
}
  \item{loc}{
%%     ~~Describe \code{loc} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{repl}{
%%     ~~Describe \code{repl} here~~
}
  \item{n.mesh}{
%%     ~~Describe \code{n.mesh} here~~
}
  \item{n.group}{
%%     ~~Describe \code{n.group} here~~
}
  \item{n.repl}{
%%     ~~Describe \code{n.repl} here~~
}
  \item{group.mesh}{
%%     ~~Describe \code{group.mesh} here~~
}
  \item{group.method}{
%%     ~~Describe \code{group.method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.spde.make.index}}
}
\examples{
loc = matrix(runif(10000*2)*1000,10000,2)
mesh = inla.mesh.create.helper(points=loc,
                               cutoff=50,
                               max.edge=c(50,500))
A = inla.spde.make.A(mesh, loc=loc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
