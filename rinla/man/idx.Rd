\name{idx}
\alias{idx}
\alias{inla.idx}

\title{Convert indexes}

\description{Convert indexes given by to triplet `(idx, group,
  replicate)' to the (one-dimensional) index used in the grouped and
  replicated model}

\usage{
inla.idx = function(
        idx,
        n = max(idx),
        group = rep(1, n),
        ngroup = max(group),
        replicate = rep(1, n),
        nrep = max(replicate))
}

\arguments{

  \item{idx}{The index within the basic model. (Legal values from `1' to `n'.)}
  \item{n}{The length `n' of the basic model.}
  \item{group}{The index within group. (Legal values from `1' to `ngroup'.)}
  \item{ngroup}{Number of groups.}
  \item{replication}{The index within replication. (Legal values from
  `1' to `nrep'.)}
  \item{nrep}{Number of replications.}
}

\value{%%
  \code{inla.idx} returns indexes in the range `1' to `n*ngroup*nrep'
  representing where the triplet `(idx,group,replicate)' is stored
  internally in the full grouped and replicated model.
}
%%

\author{Havard Rue \email{hrue@math.ntnu.no}}

\examples{
TODO
}


%  LocalWords:  idx inla ngroup nrep Havard iid
