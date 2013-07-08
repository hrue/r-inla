\name{summary.inla}
\alias{summary.inla}
\alias{summary.surv.inla}

\alias{print.summaryinla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary for a INLA fit }
\description{
Takes a fitted \code{inla} or \code{surv.inla} object produced by
\code{inla} or  \code{surv.inla} and produces
a summary from it.

}
\usage{
\method{summary}{inla}(object,...)
\method{print}{summary.inla}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{  a fitted \code{inla} object as produced by
    \code{inla}.}
  \item{x}{a \code{summary.inla} object produced by \code{summary.inla}.}
  \item{...}{ other arguments.}

}
\details{
Posterior mean and standard deviation (together with quantiles or
cdf) are printed for the fixed effects in the model.

For the random effects the function \code{summary()} prints the
posterior mean and standard deviations for the hyperparameters
}
\value{
  \code{summary.inla} returns an object of call \code{summaryinla}, a
  list with components: 
  \item{call}{the component from \code{object}.}
  \item{fixed}{the component from \code{object}.}
  \item{random}{the component from \code{object}.}
  \item{neffp}{the component from \code{object}.}          
  \item{linear.predictor}{the component from \code{object}.}
  \item{family}{the component from \code{object}.}   
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
\author{ sara martino \email{martino@math.ntnu.no}}

\seealso{ \code{\link{inla}} }
\examples{
 ## --- Continuing the Example Surg from  ?inla:

 # summary(mod.surg)

}

\keyword{ list }
