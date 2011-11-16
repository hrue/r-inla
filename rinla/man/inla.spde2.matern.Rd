\name{inla.spde2.matern}
\alias{inla.spde2.matern}
\alias{inla.spde2}
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
                  B.kappa2 = matrix(c(0,0,1),1,3),
                  sigma2.nominal = 1,
                  theta.prior.mean = NULL,
                  theta.prior.prec = NULL,
                  fractional.method = c("parsimonious", "null"))
}
\arguments{
  \item{mesh}{
    The mesh to build the model on, as an \code{\link{inla.mesh}} object.
  }
  \item{alpha}{}
  \item{B.tau}{}
  \item{B.kappa2}{}
  \item{sigma2.approx}{}
  \item{theta.prior.mean}{}
  \item{theta.prior.prec}{}
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
      \eqn{\log\kappa^2(u) = B^{\kappa^2}_0(u) + \sum_{k=1}^p B^{\kappa2}_k(u)
	\theta_k}{log kappa^2(u) = B.kappa2_0(u) + sum_{k=1}^p B.kappa2_k(u) theta_k}
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
  \code{\link{inla.mesh.generic2}}, \code{\link{inla.mesh.basis}},
  \code{\link{inla.spde.query}}
}
\examples{
n = 100
field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
loc = matrix(runif(n*2),n,2)
## One field, 2 observations per location
idx.y = rep(1:n,2)
y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))

mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
spde = inla.spde.matern(mesh)
data = list(y=y, field=mesh$idx$loc[idx.y])
formula = y ~ -1 + f(field, model=spde)
result = inla(formula, data=data, family="normal")

## Plot the mesh structure:
plot(mesh)
## For now, use old plotting methods to show fields:
col.pal = colorRampPalette(c("blue","cyan","green","yellow","red"))
plot(old.mesh.class(mesh),
     result$summary.random$field[,"mean"],
     color.palette = col.pal)
## Plot residual field:
plot(old.mesh.class(mesh),
     result$summary.random$field[,"mean"]-field.fcn(mesh$loc),
     color.palette = col.pal)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
