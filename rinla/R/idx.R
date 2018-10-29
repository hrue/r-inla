## Export: inla.idx

##!\name{idx}
##!\alias{idx}
##!\alias{inla.idx}
##!
##!\title{Convert indexes}
##!
##!\description{Convert indexes given by to triplet `(idx, group,
##!  replicate)' to the (one-dimensional) index used in the grouped and
##!  replicated model}
##!
##!\usage{
##!inla.idx(idx, n = max(idx),
##!         group = rep(1, length(idx)), ngroup = max(group),
##!         replicate = rep(1, length(idx)),  nrep = max(replicate)) 
##!}
##!
##!\arguments{
##!
##!  \item{idx}{The index within the basic model. (Legal values from `1' to `n'.)}
##!  \item{n}{The length `n' of the basic model.}
##!  \item{group}{The index within group. (Legal values from `1' to `ngroup'.)}
##!  \item{ngroup}{Number of groups.}
##!  \item{replicate}{The index within replication. (Legal values from
##!  `1' to `nrep'.)}
##!  \item{nrep}{Number of replications.}
##!}
##!
##!\value{%%
##!  \code{inla.idx} returns indexes in the range `1' to `n*ngroup*nrep'
##!  representing where the triplet `(idx,group,replicate)' is stored
##!  internally in the full grouped and replicated model.
##!}
##!%%
##!
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!
##!\examples{
##!##TODO
##!}


`inla.idx` = function(idx, n = max(idx),
                      group = rep(1, length(idx)), ngroup = max(group),
                      replicate = rep(1, length(idx)), nrep = max(replicate))
{
    ## this function might be useful to convert from (idx, group, rep)
    ## to idx, in the same way as done internally in inla.R

    stopifnot(n >= 1)
    stopifnot(ngroup >= 1)
    stopifnot(nrep >= 1)
    stopifnot(all(group >= 1))
    stopifnot(all(replicate >= 1))
    stopifnot(all(idx >= 1))
    stopifnot(all(idx <= n))
    stopifnot(ngroup >= max(group))
    stopifnot(nrep >= max(replicate))

    return (idx + (group-1)*n + (replicate-1)*n*ngroup)
}
