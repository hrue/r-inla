\name{inla.spde1.create}
\alias{inla.spde.create}
\alias{inla.spde1.create}
\alias{inla.spde1.matern}
\alias{inla.spde1.imatern}
\alias{inla.spde1.matern.osc}
\alias{inla.spde1}
\title{
  Old SPDE model objects for INLA
}
\description{
  Create an \code{inla.spde1} model object.
}
\usage{
## Create an SPDE model object:
inla.spde1.create(mesh,
                  model = c("matern", "imatern", "matern.osc"),
                  param = NULL,
                  ...)

## Shortcuts to the matern, imatern and matern.osc models:
inla.spde1.matern(mesh, ...)
inla.spde1.imatern(mesh, ...)
inla.spde1.matern.osc(mesh, ...)
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
  \item{...}{
    Additional parameters passed on to other methods.
  }
}
\details{
  Note: This is an old spde object format retained for backwards
  compatibility.  Please use \code{\link{inla.spde2}} models for new
  code.
  
  This method constructs an object for SPDE models.  Currently implemented:

  \code{model="matern"}
  
  \deqn{(\kappa^2(u)-\Delta)^{\alpha/2}(\tau(u)
    x(u))=W(u)}{(kappa^2(u)-Delta)^(alpha/2) (tau(u) x(u)) = W(u)}
  
  \code{param}:
  \itemize{
  \item{\code{alpha}}{ = 1 or 2}
  \item{\code{basis.T}}{ = Matrix of basis functions for
  \eqn{\log\tau(u)}{log tau(u)}}
  \item{\code{basis.K}}{ = Matrix of basis functions for
  \eqn{\log\kappa^2(u)}{log kappa^2(u)}}
  }

  \code{model="imatern"}
  
  \deqn{(-\Delta)^{\alpha/2}(\tau(u)
    x(u))=W(u)}{(-Delta)^(alpha/2) (tau(u) x(u)) = W(u)}
  
  \code{param}:
  \itemize{
  \item{\code{alpha}}{ = 1 or 2}
  \item{\code{basis.T}}{ = Matrix of basis functions for
  \eqn{\log\tau(u)}{log tau(u)}}
  }
}
\value{
  An \code{inla.spde1} object.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.spde2.matern}}, \code{\link{inla.mesh.2d}},
  \code{\link{inla.mesh.basis}}
}
\examples{
n = 100
field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
loc = matrix(runif(n*2),n,2)
## One field, 2 observations per location
idx.y = rep(1:n,2)
y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))

mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
spde = inla.spde.create(mesh, model="matern")
data = list(y=y, field=mesh$idx$loc[idx.y])
formula = y ~ -1 + f(field, model=spde)
result = inla(formula, data=data, family="normal")

## Plot the mesh structure:
plot(mesh)
\donttest{
if (require(rgl)) {
  ## Plot the posterior mean:
  plot(mesh, rgl=TRUE,
       result$summary.random$field[,"mean"],
       color.palette = colorRampPalette(c("blue","green","red")))
  ## Plot residual field:
  plot(mesh, rgl=TRUE,
       result$summary.random$field[,"mean"]-field.fcn(mesh$loc),
       color.palette = colorRampPalette(c("blue","green","red")))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
