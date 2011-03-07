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
          model = c("matern", "imatern"),
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
  \code{\link{inla.mesh}}, \code{\link{inla.mesh.basis}}
}
\examples{
n = 100
field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
loc = matrix(runif(n*2),n,2)
## One field, 2 observations per location
idx.y = rep(1:n,2)
y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))

mesh = inla.mesh(loc, refine=list(max.edge=0.05))
spde = inla.spde(mesh, model="matern")
data = list(y=y, field=mesh$idx$loc[idx.y])
formula = y ~ -1 + f(field, model=spde)
result = inla(formula, data=data, family="normal")

## Plot the mesh structure:
plot(mesh)
## For now, use old plotting methods to show fields:
plot(old.mesh.class(mesh),
     result$summary.random$field[,"mean"],
     color.palette = colorRampPalette(c("blue","green","red")))
## Plot residual field:
plot(old.mesh.class(mesh),
     result$summary.random$field[,"mean"]-field.fcn(mesh$loc),
     color.palette = colorRampPalette(c("blue","green","red")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
