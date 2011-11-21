\name{inla.spde2.result}
\alias{inla.spde2.result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  SPDE result extraction
}
\description{
  Exctract field and parameter values and distributions for an
  \code{inla.spde2} SPDE
  effect from an inla result object.
}
\usage{
inla.spde2.result(inla, name, spde, do.transform = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inla}{
    An \code{inla} object obtained from a call to \code{\link{inla}}
  }
  \item{name}{
    A character string with the name of the SPDE effect in the inla formula.
  }
  \item{spde}{
    The \code{inla.spde2} object used for the effect in the inla
    formula. (Note: this could have been stored in the inla output, but isn't.)
  }
  \item{do.transform}{
    If \code{TRUE}, also calculate marginals transformed to user-scale.
  }
  \item{\dots}{
    Not used.
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
mesh = inla.mesh.create.helper(points.domain=loc, max.edge=c(0.1,0.5))
spde = inla.spde2.matern(mesh)
index = inla.spde.make.index("spatial", mesh$n, n.repl=2)
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
plot(spde.result$marginals.range.nominal[[1]])
}
%% Note: the same example is used by inla.spde.make.index.Rd
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
