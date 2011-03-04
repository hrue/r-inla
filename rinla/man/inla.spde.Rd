\name{inla.spde}
\alias{inla.spde}
\alias{inla.spde.inla.mesh}
\alias{inla.spde.inla.spde}
\title{
  SPDE model objects for INLA
}
\description{
  Create and query an \code{inla.spde} model object.
}
\usage{
## Create an SPDE model object:
inla.spde(mesh,
          model = c("matern", "imatern", "heat"),
          param = NULL,
          ...)
## Query properties of an SPDE model:
inla.spde(spde, ...)
## Create an SPDE object based on INLA parameter estimates:
inla.spde(inla, name, spde, ...)
}
\arguments{
  \item{mesh}{
    The mesh to build the model on, as an \code{\link{inla.mesh}} object.
  }
  \item{model}{
    The name of the model.
  }
  \item{param}{
    Model specific parameters.
  }
  \item{spde}{
    The \code{inla.spde} object to query.
  }
  \item{...}{
    
  }
}
\details{
}
\value{
  An \code{inla.spde} object.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.mesh}}
}
\examples{
loc = matrix(runif(10*2),10,2)
y = rnorm(10)

lattice = inla.mesh.lattice(x = seq(0, 1, length.out=17),
                            y = seq(0, 1, length.out=10))
mesh = inla.mesh(loc, lattice=lattice, refine=list(max.edge=0.05))
spde = inla.spde(mesh, model="matern")
data = list(y=y, field=mesh$idx$loc)
formula = y ~ f(field, model=spde)
result = inla(formula, data=data)

plot(mesh)
plot(old.mesh.class(mesh), result$summary.random$field[,"mean"])
image(matrix(result$summary.random$field[mesh$idx$lattice,"mean"], 17, 10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
