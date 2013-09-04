\name{inla.stack}
\alias{inla.stack}
\alias{inla.stack.index}
\alias{inla.stack.data}
\alias{inla.stack.A}
\alias{inla.stack.sum}
\alias{inla.stack.join}
\alias{inla.data.stack}
\alias{inla.data.stack.info}
\alias{inla.stack.compress}
\alias{inla.stack.remove.unused}
\alias{inla.stack.LHS}
\alias{inla.stack.RHS}
\title{
  Data stacking for advanced INLA models
}
\description{
  Functions for combining data, effects and observation matrices into
  \code{inla.stack} object, and extracting information from such
  objects.
}
\usage{
## Create data stack as a sum of predictors:
inla.stack.sum(data, A, effects, tag="", compress=TRUE, remove.unused=TRUE)

## Join two or more data stacks:
inla.stack.join(..., compress=TRUE, remove.unused=TRUE)

## Shorthand for inla.stack.join and inla.stack.sum:
inla.stack(..., compress=TRUE, remove.unused=TRUE)

## Compress an existing stack:
inla.stack.compress(stack, remove.unused=TRUE)

## Remove unused entries from an existing stack:
inla.stack.remove.unused(stack)

## Extract tagged indices:
inla.stack.index(stack, tag)

## Extract data for inla call, and optionally join with other variables:
inla.stack.data(stack, ...)

## Extract "A matrix" for control.predictor:
inla.stack.A(stack)

## Extract data associated with the "left hand side" of the model:
## (e.g. the data itself, Ntrials, link, E, ...)
inla.stack.LHS(stack)

## Extract data associated with the "right hand side" of the model
## (all the covariates/predictors)
inla.stack.RHS(stack)
}
\arguments{
  \item{data}{
    A list of data vectors.
  }
  \item{A}{
    A list of observation matrices.
  }
  \item{effects}{
    A collection of effects/predictors.  Each list element corresponds
    to an observation matrix, and must either be a single vector or a
    list of vectors.
  }
  \item{tag}{
    A string specifying a tag for later identification.
  }
  \item{compress}{
    If \code{TRUE}, compress the model by removing duplicated rows of
    effects, replacing the corresponding A-matrix columns with a single
    column containing the sum.
  }
  \item{remove.unused}{
    If \code{TRUE}, compress the model by removing rows of effects
    corresponding to all-zero columns in the A matrix (and removing
    those columns).
  }
  \item{stack}{
    An \code{inla.data.stack} object, created by a call to \code{inla.stack},
    \code{inla.stack.sum}, or \code{inla.stack.join}.
  }
  \item{...}{
    For \code{inla.stack.join}, two or more data stacks of class
    \code{inla.data.stack}, created by a call to \code{inla.stack},
    \code{inla.stack.sum}, or \code{inla.stack.join}.

    For \code{inla.stack.data}, a list of variables to be joined with
    the data list.
  }
}
\details{
  \itemize{
    \item[\eqn{n_l}{n_l}] effect blocks
    \item[\eqn{n_k}{n_k}] effects
    \item[\eqn{n_i}{n_i}] data values
    \item[\eqn{n_{j,l}}{n_jl}] effect size for block \eqn{l}{l}
    \item[\eqn{n_j}{n_j}] \eqn{= \sum_{l=1}^{n_l} n_{j,l}}{sum_l n_jl} total effect size
  }
  Input:
  \itemize{
    \item[\code{data}] \eqn{(y^1, \ldots, y^p)}{} \eqn{p}{p} vectors,
    each of length \eqn{n_i}{n_i}
    \item[\code{A}] \eqn{(A^1, \ldots, A^{n_l})}{} matrices of size
    \eqn{n_i \times n_{j,l}}{n_i by n_jl}
    \item[\code{effects}] \eqn{\left((x^{1,1},\ldots,x^{n_k,1}), \ldots,
      (x^{1,n_l},\ldots,x^{n_k,n_l})\right)}{} collections of effect vectors
    of length \eqn{n_{j,l}}{n_jl}
  }
  
  \deqn{
    \mbox{predictor}(y^1, \ldots, y^p) \sim
    \sum_{l=1}^{n_l} A^l \sum_{k=1}^{n_k} g(k, x^{k,l})
    = \tilde{A} \sum_{k=1}^{n_k} g(k, \tilde{x}^k)
  }{
    predictor(y^1, \ldots, y^p)
    ~ sum_{l=1}^{n_l} A^l sum_{k=1}^{n_k} g(k, x^{k,l})
    = tilde{A} sum_{k=1}^{n_k} g(k, tilde{x}^k)
  }
  where
  \deqn{
    \tilde{A} = \mbox{cbind}\left( A^1, \ldots, A^{n_l} \right)
  }{
    tilde{A} = cbind( A^1, ..., A^{n_l} )
  }
  \deqn{
    \tilde{x}^k = \mbox{rbind}\left( x^{k,1}, \ldots, x^{k,n_l} \right)
  }{
    tilde{x}^k = rbind( x^{k,1}, ..., x^{k,n_l} )
  }
  and for each block \eqn{l}{l}, any missing \eqn{x^{k,l}} is replaced by an
  \code{NA} vector.
}
\value{A data stack of class \code{inla.data.stack}.  Elements:
  \itemize{
    \item{\code{data} }{\eqn{=(y^1, \ldots, y^p, \tilde{x}^1, \ldots, \tilde{x}^{n_k})}}
    \item{\code{A} }{\eqn{=\tilde{A}}}
    \item{\code{data.names} }{List of data names, length \eqn{p}}
    \item{\code{effect.names} }{List of effect names, length \eqn{n_k}}
    \item{\code{n.data} }{Data length, \eqn{n_i}}
    \item{\code{index} }{List indexed by \code{tag}s, each element indexing
    into \eqn{i=1, \ldots, n_i}}
  }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.spde.make.A}},
  \code{\link{inla.spde.make.index}}
}
\examples{
n = 200
loc = matrix(runif(n*2),n,2)
mesh = inla.mesh.create.helper(points.domain=loc,
                               max.edge=c(0.05, 0.2))
proj.obs = inla.mesh.projector(mesh, loc=loc)
proj.pred = inla.mesh.projector(mesh, loc=mesh$loc)
spde = inla.spde2.matern(mesh,
                         B.tau=cbind(log(1), 1, 0),
                         B.kappa=matrix(c(log(sqrt(8)/0.2), 0, 1), 1, 3))

covar = rnorm(n)
field = inla.qsample(n=1, Q=inla.spde2.precision(spde, theta=c(0,0)))[,1]
y = 2*covar + inla.mesh.project(proj.obs, field)

A.obs = inla.spde.make.A(mesh, loc=loc)
A.pred = inla.spde.make.A(mesh, loc=proj.pred$loc)
stack.obs =
    inla.stack(data=list(y=y),
               A=list(A.obs, 1),
               effects=list(c(list(intercept=rep(1, mesh$n)),
                              inla.spde.make.index("spatial", spde$n.spde)),
                            covar=covar),
               tag="obs")
stack.pred =
    inla.stack(data=list(y=NA),
               A=list(A.pred),
               effects=list(c(list(intercept=rep(1, mesh$n)),
                              inla.spde.make.index("spatial", mesh$n))),
               tag="pred")
stack = inla.stack(stack.obs, stack.pred)

formula = y ~ -1 + intercept + covar + f(spatial, model=spde)
result1 = inla(formula,
               data=inla.stack.data(stack.obs, spde=spde),
               family="gaussian",
               control.predictor=list(A=inla.stack.A(stack.obs), compute=TRUE),
               verbose=TRUE)

plot(y,result1$summary.fitted.values[inla.stack.index(stack.obs,"obs")$data,"mean"])

result2 = inla(formula,
               data=inla.stack.data(stack, spde=spde),
               family="gaussian",
               control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
               verbose=TRUE)

field.pred = inla.mesh.project(proj.pred,
    result2$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"])

dev.new()
plot(field, field.pred)
dev.new()
image(inla.mesh.project(mesh,
    field=field,
    dims=c(200,200)))
dev.new()
image(inla.mesh.project(mesh,
    field=field.pred,
    dims=c(200,200)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
