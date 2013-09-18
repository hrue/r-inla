\name{inla.spde.make.index}
\alias{inla.spde.make.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  SPDE model index vector generation
}
\description{
  Generates a list of named index vectors for an SPDE model.
}
\usage{
inla.spde.make.index(name,
                     n.spde,
                     n.group = 1,
                     n.repl = 1,
                     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
    A character string with the base name of the effect.
}
\item{n.spde}{
  The size of the model, typically from \code{spde$n.spde}.
}
\item{n.group}{
  The size of the \code{group} model.
}
\item{n.repl}{
  The number of model replicates.
}
\item{...}{
  Additional parameters.  Currently unused.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of named index vectors.
\item{name }{Indices into the vector of latent variables}
\item{name.group }{'group' indices}
\item{name.repl }{Indices for replicates}
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
\code{\link{inla.spde.make.A}},
\code{\link{inla.spde2.result}}
}
\examples{
loc = matrix(runif(100*2),100,2)
mesh = inla.mesh.create.helper(points.domain=loc, max.edge=c(0.1,0.5))
spde = inla.spde2.matern(mesh)
index = inla.spde.make.index("spatial", spde$n.spde, n.repl=2)
spatial.A = inla.spde.make.A(mesh, loc,
                             index=rep(1:nrow(loc), 2),
                             repl=rep(1:2, each=nrow(loc)))
y = 10+rnorm(100*2)
stack = inla.stack(data=list(y=y),
                   A=list(spatial.A),
                   effects=list(c(index, list(intercept=1))),
                   tag="tag")
data = inla.stack.data(stack, spde=spde)
formula = y ~ -1 + intercept + f(spatial, model=spde,
                                 replicate=spatial.repl)
result = inla(formula, family="gaussian", data=data, 
              control.predictor=list(A=inla.stack.A(stack)))
spde.result = inla.spde2.result(result, "spatial", spde)
}
%% Note: the same example is used by inla.spde2.result.Rd
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
