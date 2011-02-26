\name{inla.models.general}
\alias{inla.models.general}

\title{Specifying generalized additive models to be solved using
  INLA}
\description{
  This page is intended to provide some more information on how to
  specify generalized additive (GAM) models that should be solved using
  \code{\link{inla}}.
  Assume that, in general, we are interested in
  modelling some response variable \eqn{y_i}{y} with an exponential
  family distribution (e.g. Normal, Poisson, binomial, etc.) using predictor variables
  \eqn{x_{1i}, x_{2i}, \ldots}{x1,x2,...}, and let \eqn{\mu_i \equiv
    E(y_i)}{m=E(y)}. A typical model might be:
  \deqn{g(\mu_i) = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} +
    w_{1i}\ f_1(x_{3i})+w_{2i}\ f_2(x_{4i})+\mbox{offset}_i}{g(m)=b0 +
  b1x1 + b2x2 + w1 f1(x3) +
    w2 f2(x4)+offset }
  where \eqn{g} is a smooth monotonic `link' function,  \eqn{f_1}{f1}
  and \eqn{f_2}{f2} are smooth functions of the respective covariate
  values, and \eqn{\mbox{offset}_i}{offset}
  is can be used to specify an a priori known component to
  be included in the linear predictor during fitting.
The smooth functions of covariates \eqn{f_1}{f1}
  and \eqn{f_2}{f2} might be weigthed through the know weight vectors
  \eqn{\mathbf{w}_1}{w1} and \eqn{\mathbf{w}_2}{w2}. 
    
  The parameters \eqn{\beta_1}{b1},\eqn{\beta_2}{b2} etc. express the
  effect on  \eqn{g(\mu)}{g(m)} of the covariates \eqn{x_1}{x1}, \eqn{x_2}{x2}
  assumed to have a \emph{fixed effect}.
  The functions \eqn{f_1}{f1},\eqn{f_2}{f2} etc. can be used to represent
  smooth functions of covarites (as in GAM models), smooth effect of time
  (as in time series models) or even spatial effect (as in the so called
  Geoadiitive effects). 
  The key idea here is that the
  dependence of the response on the predictors can be represented as a
  parametric sub-model plus the sum of some smooth functions of one
  predictor variables.
  Thus the model is quite flexible
  relative to strictly parametric linear or generalized linear models.

  We assume a Bayesian perspective therefore all unknowm quantities in
  the model have a prior distribution. The parameters for the fixed
  effect covariates \eqn{\beta_1}{b1},\eqn{\beta_2}{b2} are assumed to
  have a centered Gaussian distribution with \emph{known} precision parameter

  \deqn{\beta_i\sim N(0,\mbox{prec}^{-1})}{b\sim N(0,prec^(-1))}

  As for the smooth functions \eqn{f_k}{fk} are modelled by assuming that
  the corresponding vector of function evaluations
  \eqn{\mathbf{f}_k=\{f_{1k},\ldots,f(n_k
    k)\}}{\bold{f}_k=(f_1k,\ldots,f_nk)} has a cenered Gaussian
  distribution with precision matrix governed by some hyperparameters \eqn{\theta}{\theta}

  \deqn{\mathbf{f}_k\sim N(\mathbf{0},\mathbf{Q}(\theta)}{\bold{f}_k\sim
    N(0,Q(\theta))}

  This Gaussian model is often intrinsic.

  A last level on the model consists in the specification of prior
  distribution for the hyperparameters of the model, namely \eqn{\theta}{\theta}
  and, possibly, some addtional hyperparameters in the likelihood.
  
  Note one important point. In order for the model to be
  identifiable, if the model for  the smooth functions is intrinsic, \emph{sum
  to zero constraint} have to be added. Such
  constraints are  applied as default by \code{inla}.

  Specification of the distribution and link function is done using
  the \code{family} argument to \code{inla} and works in the same way as for
  \code{glm}. This page hence assumes the default identity link normal
  error family, since the generalizations are easy.

  On the other side, specifications of the prior model for the smooth
  functions is done using the function \code{\link{f}}, for example
  the formula:\cr

  \code{y~ x1+x2+f(x3,model="iid")+f(x4,model="iid")}

  specifies both \eqn{x_3}{x3} and \eqn{x_4}{x4} to be random effects,
  while the formula

  \code{y~ x1+x2+f(x3,model="rw2")+f(x4,model="iid")}

  specifies \eqn{x_3}{x3} to have a random walk of order 2 prior and \eqn{x_4}{x4} to be random effect.


  
}
\author{sara martino \email{martino@math.ntnu.no}}
\seealso{\code{\link{inla}} ,  \code{\link{f}} }

\keyword{models}
\keyword{regression}

 
