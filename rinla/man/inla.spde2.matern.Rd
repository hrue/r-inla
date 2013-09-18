\name{inla.spde2.matern}
\alias{inla.spde2.matern}
\title{
  Matern SPDE model object for INLA
}
\description{
  Create an \code{inla.spde2} model object for a Matern model.
}
\usage{
inla.spde2.matern(mesh,
                  alpha = 2,
                  B.tau = matrix(c(0,1,0),1,3),
                  B.kappa = matrix(c(0,0,1),1,3),
                  prior.variance.nominal = 1,
                  prior.range.nominal = NULL,
                  prior.tau = NULL,
                  prior.kappa = NULL,
                  theta.prior.mean = NULL,
                  theta.prior.prec = NULL,
                  fractional.method = c("parsimonious", "null"))
}
\arguments{
  \item{mesh}{
    The mesh to build the model on, as an \code{\link{inla.mesh}} object.
  }
  \item{alpha}{Fractional operator order, \eqn{0<\alpha\leq 2} supported. (\eqn{\nu=\alpha-d/2})}
  \item{B.tau}{Matrix with specification of log-linear model for \eqn{\tau}.}
  \item{B.kappa}{Matrix with specification of log-linear model for \eqn{\kappa}.}
  \item{prior.variance.nominal}{Nominal prior mean for the field variance}
  \item{prior.range.nominal}{Nominal prior mean for the spatial range}
  \item{prior.tau}{Prior mean for tau (overrides \code{prior.variance.nominal})}
  \item{prior.kappa}{Prior mean for kappa (overrides \code{prior.range.nominal})}
  \item{theta.prior.mean}{(overrides \code{prior.*})}
  \item{theta.prior.prec}{Scalar, vector or matrix, specifying the joint prior
  precision for \eqn{theta}.}
  \item{fractional.method}{}
}
\details{
  This method constructs a Matern SPDE model.

  \deqn{(\kappa^2(u)-\Delta)^{\alpha/2}(\tau(u)
    x(u))=W(u)}{(kappa^2(u)-Delta)^(alpha/2) (tau(u) x(u)) = W(u)}
  
  \itemize{
    \item{\eqn{\alpha}{alpha}}{ = 2 (code for other values is in progress)}
    \item{}{
      \eqn{\log\tau(u) = B^\tau_0(u) + \sum_{k=1}^p B^\tau_k(u)
	\theta_k}{log tau(u) = B.tau_0(u) + sum_{k=1}^p B.tau_k(u) theta_k}
    }
    \item{}{
      \eqn{\log\kappa(u) = B^{\kappa}_0(u) + \sum_{k=1}^p B^{\kappa}_k(u)
	\theta_k}{log kappa(u) = B.kappa_0(u) + sum_{k=1}^p B.kappa_k(u) theta_k}
    }
  }
}
\value{
  An \code{inla.spde2} object.
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
  \code{\link{inla.mesh.2d}}, \code{\link{inla.mesh.create}},
  \code{\link{inla.mesh.1d}}, \code{\link{inla.mesh.basis}},
  \code{\link{inla.spde2.generic}}
}
\examples{
n = 100
field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
loc = matrix(runif(n*2),n,2)
## One field, 2 observations per location
idx.y = rep(1:n,2)
y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))

mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
spde = inla.spde2.matern(mesh)
data = list(y=y, field=mesh$idx$loc[idx.y])
formula = y ~ -1 + f(field, model=spde)
result = inla(formula, data=data, family="normal")

## Plot the mesh structure:
plot(mesh)
\donttest{
if (require(rgl)) {
  col.pal = colorRampPalette(c("blue","cyan","green","yellow","red"))
  ## Plot the posterior mean:
  plot(mesh, rgl=TRUE,
       result$summary.random$field[,"mean"],
       color.palette = col.pal)
  ## Plot residual field:
  plot(mesh, rgl=TRUE,
       result$summary.random$field[,"mean"]-field.fcn(mesh$loc),
       color.palette = col.pal)
}
}
result.field = inla.spde.result(result, "field", spde)
plot(result.field$marginals.range.nominal[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
