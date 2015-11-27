## Export: inla.scale.model

##! \name{scale.model}
##! \alias{inla.scale.model}
##! \alias{scale.model}
##!
##! \title{Scale an intrinsic GMRF model}
##!
##! \description{This function scales an intrinsic GMRF model so the geometric mean of the
##!              marginal variances is one}
##!
##! \usage{
##!     inla.scale.model(Q, constr, eps = sqrt(.Machine$double.eps))
##! }
##!
##! \arguments{
##!   \item{Q}{A SPD matrix,  either as a (dense) matrix or \code{sparseMatrix}}
##!   \item{constr}{Linear constraints spanning the null-space of \code{Q};
##!                 see \code{?INLA::f} and argument \code{extraconstr}} 
##!   \item{eps}{A small constant added to the diagonal of \code{Q} to name it non-singular}
##!  }
##! \value{
##!   \code{inla.scale.model} returns a \code{sparseMatrix} of type \code{dgTMatrix} 
##!   scaled so the geometric mean of the marginal variances (of the non-singular part of
##!   \code{Q}) is one.
##! }
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
##!
##! \examples{
##! ## make an intrinsic GMRF (model="besag")
##! data(Germany)
##! g = system.file("demodata/germany.graph", package="INLA")
##! Q = -inla.graph2matrix(g)
##! diag(Q) = 0
##! diag(Q) = -rowSums(Q)
##! stopifnot(all(rowSums(Q) == 0))
##! n = dim(Q)[1]
##! Q.scaled = inla.scale.model(Q, constr = list(A = matrix(1, 1, n), e=0))
##! }

inla.scale.model = function(Q, constr, eps = sqrt(.Machine$double.eps))
{
    if (missing(constr)) {
        stop("Argument missing: constr = list(A=matrix(...), e=c(...))")
    }
    Q = inla.as.sparse(Q)
    n = dim(Q)[1]
    res = inla.qinv(Q + Diagonal(n) * max(diag(Q)) * eps, constr = constr)
    fac = exp(mean(log(diag(res))))
    Q = fac * Q

    return (Q)
}
