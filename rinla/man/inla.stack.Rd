\name{inla.stack}
\alias{inla.stack}
\alias{inla.stack.default}
\alias{inla.stack.inla.data.stack}
\alias{inla.data.stack}
\title{
  Data stacking for advanced INLA models
}
\description{
  Combine data, effects and observation matrices
}
\usage{
## Create data stack:
inla.stack(data, A, effects, tag = NULL, strict = TRUE, ...)

## Join two or more data stacks:
inla.stack(...)
}
\arguments{
  \item{data}{
    A list of data vectors.
  }
  \item{A}{
    A list of observation matrices.
  }
  \item{effects}{
    A list of effects.  Each list element corresponds to an observation
    matrix, and must either be a single vector or a list of vectors.
  }
  \item{tag}{
    A string specifying a tag for later identification.
  }
  \item{strict}{
    If \code{TRUE}, all input variables must have matching sizes (length of
    data = rows of A matrices, length of effects = columns of A matrices).
    If \code{FALSE}, an attempt is made to pad short variables with \code{NA}.
  }
  \item{...}{
    Two or more data stacks of class \code{inla.data.stack}, created by
    a call to \code{inla.stack}.
    (Used only when joining data stacks, not when creating them.)
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
  \code{\link{inla.spde.make.index}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ fmesher }
