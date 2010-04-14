\name{f}
\alias{f}
\title{ Defining smooth and spatial effects in INLA formulae }
\description{

  Function used in definition of smooth and spatial terms within \code{inla} model
  formulae. The function does not evaluate a anything - it
  exists purely to help set up a model. The function specifies one
  smooth function in the linear predictor (see \code{\link{inla.models}}) as
  \deqn{w\ f(x)}{weight*f(var)}
  
}
\usage{
`f` =
function(...,
         model=NULL,
         n=NULL,
         nrep=NULL,
         replicate=NULL,
         initial=NULL,
         season.length=NULL,
         prior=NULL,
         param = NULL,
         fixed = NULL,
         constr = NULL,
         extraconstr=NULL,
         values=NULL,
         cyclic = NULL,
         diagonal = NULL, 
         graph.file=NULL,
         cdf=NULL,
         quantiles=NULL,
         Cmatrix=NULL,
         rankdef=NULL,
         Z = NULL,
         nrow = NULL,
         ncol = NULL,
         nu = NULL,
         mean.linear=NULL,
         prec.linear=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{...}{Name of the covariate and, possibly of the weights
  vector. NB: order counts!!!! The first specified term is the covariate
  and the second one is the vector of weights.}
  
  \item{model}{Model for the smooth or spatial effect. Available
    choices are:\cr
    \tabular{lll}{
      
      Name \tab Model \tab N. of hyperparameters \cr
      \code{"iid"} \tab random effect \tab 1\cr
      \code{"rw1"} \tab Random walk of order 1 \tab 1\cr
      \code{"rw2"} \tab Random walk of order 2 \tab 1\cr
      \code{"besag"} \tab Spatial conditional autoregressive \tab 1\cr
      \code{"ar1"} \tab Autoregressice of order 1 \tab 2\cr
      \code{"seasonal"} \tab Model for Seasonal component \tab 1\cr
      \code{"generic0"} \tab Generic model \tab 1\cr%
      \code{"generic1"} \tab Generic model \tab 1\cr%
      \code{"linear"}\tab special model for linear effect with different
      priors for the mean and precision \tab 2\cr
      \code{"2diid"}\tab 2 dimensional iid (see inla manual for details)\tab 3\cr
      \code{"2diidwishart"}\tab 2 dimensional iid with Wishart-prior (see inla manual for details)\tab 3\cr
      \code{"3diidwishart"}\tab 3 dimensional iid with Wishart-prior (see inla manual for details)\tab 6\cr
      \code{"rw2d"}\tab Two-dimentional 2.order RW-model\tab 1\cr
      \code{"matern2d"}\tab Two-dimentional Matern-fields\tab 2\cr
    }
  }

  \item{n}{An optional argument which defines the dimension of the model
    if this is different from \code{length(sort(unique(covariate)))}}
    
  \item{nrep}{(DO NOT USE) An optional argument which defines how many
    times this model should replicate itself. This is an expert option;
    see examples (which are still missing...) for usage.}

  \item{replicate}{(DO NOT USE) An optional argument which defines
    which replicate the covariate value is intended for. This is an
    expert option; see examples (which are still missing...) for usage.}
    
  \item{initial}{Vector of length 1 or 2 indicating the starting values
    for the optimization algorithm. If \code{fixed=T} the value ti which
    the parameters are fixed is determines through \code{initial}. (See Details) }
  
  \item{season.length}{Lenght of the seasonal compoment (ONLY if
    \code{model="seasonal"})  }

  \item{prior}{Prior distribution(s) for the hyperparameters of the
    random model. The default value depends on the type of model.}

  \item{param}{Vector of lenght 2 or 4 indicating the parameters \eqn{a}{a} and \eqn{b}{b}
    of the prior distribution for the hyperparameters (See Details) }
  
  \item{fixed}{Vector of length 1 or 2 indicating wheater the
    hyperparameters of the model are fixed or random.  (See Details)}

  \item{constr}{A boolean variable indicating whater to set a sum to 0
    constraint on the term. By default the sum to 0 constraint is
    imposed on all intrinsic models ("iid","rw1","rw1","besag", etc..).}

  \item{extraconstr}{This argument defines extra linear constraints. The
    argument is a list with two elements, a matrix \code{A} and a vector
    \code{e}, which defines the extra constraint \code{Ax = e}; for
    example \code{extraconstr = list(A = A, e=e)}. The number of columns
    of \code{A} must correspond to the length of this \code{f}-model.
    Note that this constraint comes additional to the sum-to-zero
    constraint defined if \code{constr = TRUE}.}

  \item{values}{An optional vector giving all values assumed by the
  covariate for which we want estimated the effect. It must be a numeric
  vector, a vector of factors or \code{NULL}.}
  
  \item{cyclic}{A boolean specifying wheather the model is
    cyclical. Only valid for "rw1" and "rw2" models, is cyclic=T then
    the sum to 0 constraint is removed. For the correct form of the grah
    file see \cite{Martinoand Rue (2008)}.}

  \item{graph.file}{Name of the file containing the graph of the model.
    For the correct for of the graph file see \cite{Martino and Rue (2008)} }

  \item{cdf}{A vector of maximum 10 values between 0 and 1
    \eqn{x(0),x(1),\ldots}{x(0),x(1),\ldots}. The function returns, for
    each posterior marginal the probabilities
    \deqn{\mbox{Prob}(X<x(p))}{Prob(X<x(p))} }

  \item{quantiles}{A vector of maximum 10 quantiles,
    \eqn{p(0),p(1),\dots}{p(0),p(1),\ldots} to compute for each
    posterior marginal. The function returns, for each posterior
    marginal, the values \eqn{x(0),x(1),\dots}{x(0),x(1),\ldots} such
    that \deqn{\mbox{Prob}(X<x(p))=p}{Prob(X<x)=p} }

  \item{Cmatrix}{The specification of the precision matrix for the
    generic models (up to a scaling constant), and is only used if
    \code{model="generic0"} or \code{model="generic1"}.  \code{Cmatrix}
    is a list of type \code{Cmatrix = list(i = c(), j = c(), Cij =
    c())}, where \code{i}, \code{j} and \code{Cij} are vectors of the
    non-zero elements of \code{C}. Note that \code{i} starts \code{1},
    and only the upper or lower part of \code{C} has to be given.

    Alternatively, \code{Qmatrix} is the name of a file giving the
    precision matrix as described in the \code{inla}-manual. }

  \item{rankdef}{A number \bold{defining} the rank deficiency of the
    model, with sum-to-zero constraint and possible extra-constraints
    taken into account. See details.  }

  \item{nrow}{Number of rows for 2d-models}

  \item{ncol}{Number of columns for 2d-models}

  \item{nu}{Smoothing parameter for the Matern2d-model, nu in \code{c(0,1,2,3)}}

  \item{mean.linear}{Prior mean for the linear component, only used if
    \code{model="linear"}}
  
  \item{prec.linear}{Prior precision for the linear component, only used
    if \code{model="linear"}}
  
}
\details{
   All available models, except for "ar1", contain only one
  hyperparamter, namely the log-precision parameter
  \eqn{\lambda}{\lambda}. The prior distribution for the log-precision is
  LogGamma with parameters  \eqn{a}{a} and \eqn{b}{b}. The LogGamma
  distribution is defined such that,
  \eqn{X\sim\mbox{LogGamma}(a,b)}{X~LogGamma(a,b)} iff
  \eqn{Y=\exp(X)\sim \mbox{Gamma}(a,b)}{Y=exp(X)~Gamma(a,b)} with
  mean \eqn{a/b}{a/b} and variance \eqn{a/b^2}{a/b^2}.\cr
  For all above-mentioned models the argument \code{param} takes the
  form \code{c(a,b)}, argument \code{initial} indicates the starting
  value for the log-precision in the optimization algorithm and the
  argument \code{fixed} indicates wheather or not the precision
  \eqn{\lambda}{\lambda} is fixed. The  argument \code{initial} contains
  one scalar number which is the starting value for log-precision.
  
  The "ar1" model is defined as
  \deqn{x_t=\phi x_{t-1}+\epsilon_t\ t=1,\ldots,n}{x_t=phi*x_{t-1}+\epsilon_t\ t=1,\ldots,n}
  where
  \eqn{\epsilon_t\sim\mathcal{N}(0,\lambda^{-1/2})}{\epsilon_t~N(0,1/\lambda)}, and has therefore two hyperparameters
  \eqn{(\lambda,\phi)}{(lambda,phi)}.\cr
  The prior for the log-precision parameter \eqn{\lambda}{lambda} is a
  \eqn{\mbox{LogGamma}(a,b)}{LogGamma(a,b)} as before.
  As for the hyperparameter \eqn{\phi}{phi}, we consider the
  transformation
  \deqn{\kappa=\mbox{logit}\frac{\phi+1}{2}}{kappa=logit((phi+1)/2)}
  and assign a 
  a 0 mean
  normal with precision parameter \eqn{\mbox{prec}_k}{prec_k} as prior
  to \eqn{\kappa}{kappa}. \cr
  For the "ar1" model the argument \code{param} takes the
  form \code{c(a,b,0,prec_k)}, argument \code{fixed} is a vector of
  lenght 2, the first element indicates wheather \eqn{\lambda}{lambda}
  is fixed and the second element wheather \eqn{\phi}{phi} is fixed. 
  The  argument \code{initial} contains two numbers, the first is the
  starting value for the log-precision and the second is the starting
  value for the parameters \eqn{\kappa}{kappa} defined above.

  
  There is no default value for \code{rankdef}, if it is not defined by
  the user  then it is computed by the
  rank deficiency of the prior model (for the generic model, the
  default is zero), plus 1 for the sum-to-zero constraint if the
  prior model is proper, plus the number of extra
  constraints. \bold{Oops:} This can be wrong, and then the user
  must define the \code{rankdef} explicitely.
}
\value{
  The function return a list of objects.
  
  (NOTE: the description of the output is not accurate since it is used
  only for internal computations.)
}
\references{
  Rue, H. and Martino, S. and Chopin, N. (2008)
  \emph{Approximate Bayesian Inference for latent Gaussian models using
    Integrated Nested Laplace Approximations}

  Rue, H and Held, L. (2005)
  \emph{Gaussian Markov Random Fields - Theory and Applications}
  Chapman and Hall
  
  Martino, S. and Rue, H. (2008)
  \emph{Implementing Approximate Bayesian Inference using Integrated
    Nested Laplace Approximation: a manual for the inla program}
  Preprint N.2 from Dep. of Mathematical Sciences (NTNU Norway)
}
\author{ sara martino \email{martino@math.ntnu.no} }

\seealso{\code{\link{inla.models}}, \code{\link{inla}},  \code{\link{hyperpar.inla}}}


\keyword{ models }

