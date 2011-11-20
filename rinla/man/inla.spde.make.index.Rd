\name{inla.spde.make.index}
\alias{inla.spde.make.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inla.spde.make.index(name,
                     n.mesh,
                     n.group = 1,
                     n.repl = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
%%     ~~Describe \code{name} here~~
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
  \item{n.field}{
%%     ~~Describe \code{n.field} here~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
loc = matrix(runif(100*2),100,2)
mesh = inla.mesh.create.helper(points.domain=loc, max.edge=c(0.05,0.5))
index = inla.spde.make.index("spatial", mesh$n, n.repl=2)
spatial.A = inla.spde.make.A(mesh, loc,
                             index=rep(1:nrow(loc), 2),
                             repl=rep(1:2, each=nrow(loc)))
y = 10+rnorm(100*2)
stack = inla.stack(data=list(y=y),
                   A=spatial.A,
                   effects=c(index, list(intercept=1)),tag="tag")
data = inla.stack.data(stack)
formula = y ~ -1 + intercept + f(spatial, model=spde,
                                 replicate=spatial.repl)
result = inla(formula, family="gaussian", data=data, 
              control.predictor(A=inla.stack.A(stack)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
