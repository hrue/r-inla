## Export: inla.knmodels.sample

##! \name{inla.knmodels.sample}
##! \alias{inla.knmodels.sample}
##! \title{Spacetime interaction models sampler function}
##! \description{
##!    It implements the sampling method for the models 
##!    in Knorr-Held, L. (2000) considering the algorithm 
##!    3.1 in Rue & Held (2005) book.
##! } 
##! \usage{
##!inla.knmodels.sample(
##!  graph,
##!  m,
##!  type=4,
##!  intercept=0, 
##!  tau.t=1,
##!  phi.t=0.7,
##!  tau.s=1,
##!  phi.s=0.7,
##!  tau.st=1,
##!  ev.t=NULL,
##!  ev.s=NULL)
##!}
##! \arguments{
`inla.knmodels.sample` =
    function(
        ##! \item{graph}{}
        graph,
        ##! \item{m}{Time dimention.}
        m,
        ##! \item{type}{Integer from 1 to 4 to identify one
        ##!    of the four interaction type.}
        type=4,
        ##! \item{intercept}{A constant to be added to the linear predictor}
        intercept=0, 
        ##! \item{tau.t}{Precision parameter for the main temporal effect.}
        tau.t=1,
        ##! \item{phi.t}{Mixing parameter in the \code{bym2} model 
        ##!    assumed for the main temporal effect.}
        phi.t=0.7,
        ##! \item{tau.s}{Precision parameter for the main spatial effect.}
        tau.s=1,
        ##! \item{phi.s}{Mixing parameter in the \code{bym2} model
        ##!    assumed for the main spatial effect.}
        phi.s=0.7,
        ##! \item{tau.st}{Precision parameter for the spacetime effect.}
        tau.st=1,
        ##! \item{ev.t}{Eigenvalues and eigenvectors of the temporal
        ##!    precision matrix structure.}
        ev.t=NULL,
        ##! \item{ev.s}{Eigenvalues and eigenvectors of the spatial 
        ##!    precision matrix structure.}
        ev.s=NULL)    
##! }
##! \value{
##!  A list with the following elements
##!   \item{time}{The time index for each obervation, 
##!       with length equals m*n.} 
##!   \item{space}{The spatial index for each obervation, 
##!       with length equals m*n.} 
##!   \item{spacetime}{The spacetime index for each obervation, 
##!       with length equals m*n.} 
##!   \item{x}{A list with the following elements}
##!   \item{t.iid}{The unstructured main temporal effect part.}
##!   \item{t.str}{The structured main temporal effect part.}
##!   \item{t}{The main temporal effect with length equals 2m.}
##!   \item{s.iid}{The unstructured main spatial effect part.}
##!   \item{s.str}{The structured main spatial effect part.}
##!   \item{s}{The main spatial effect with length equals 2n.}
##!   \item{st}{The spacetime interaction effect with length equals m*n.}
##!   \item{eta}{The linear predictor with length equals n*m.}
##! }
##! \author{Elias T. Krainski}
##! \seealso{
##!     \code{\link{inla.knmodels}} for model fitting
##! }
{
    type <- pmatch(type, 1:4)
    if (!any(type==(1:4))) stop("'type' must be 1, 2, 3 or 4!") 
    qsample <- function(q, ev=NULL, verbose=FALSE) {
        ## algorithm 3.1 to sample from a precision matrix 
        ## or from the eigenvalue/vector pairs when suplied 
        if (is.null(ev)) 
            ev <- eigen(as.matrix(q))
        n <- length(ev$values)
        jj <- which(ev$values>sqrt(.Machine$double.eps))
        k <- n-length(jj)
        if (verbose) cat('rankdef =', k, '\n')
        y <- rnorm(n-k, 0, 1/sqrt(ev$values[jj]))
        return(colSums(t(ev$vectors[, jj])*y))
    }
    if (missingArg(ev.t)) {
        ev.t <- eigen(as.matrix(crossprod(diff(
            Diagonal(m), lag=1, differences=1))))
    } else m <- length(ev.t$values)
    if (missingArg(ev.s)) {
        if (missingArg(graph)) 
            stop("'graph' or 'ev.s' must be provided!")        
        graph <- inla.graph2matrix(graph) 
        n <- nrow(graph)
        R.s <- inla.scale.model(Diagonal(n, colSums(graph)) - graph,
            constr=list(A=matrix(1, 1, n), e=0))
        ev.s <- eigen(as.matrix(R.s))
    } else n <- length(ev.s$values)
    dat <- list(time=rep(1:m, each=n), space=rep(1:n, m), spacetime=1:(m*n), x=list()) 
    ### AR(1), as used before:
    ##    dat$x$t <- arima.sim(model=list(ar=rho), n=m, ### sample with marginal variance = 1/tau.t
    ##                         rand.gen=function(leng) rnorm(leng, 0, sqrt((1-rho^2)/tau.t)))
    dat$x$t.str <- qsample(ev=ev.t)
    dat$x$t.iid <- rnorm(m, 0.0, 1.0)
    dat$x$time <- c((sqrt(phi.t)*dat$x$t.str + sqrt(1-phi.t)*dat$x$t.iid)/sqrt(tau.t),
                    dat$x$t.str)
    dat$x$s.iid <- rnorm(length(ev.s$value), 0, 1) 
    dat$x$s.str <- qsample(ev=ev.s)
    dat$x$space <- c((sqrt(phi.s)*dat$x$s.str + sqrt(1-phi.s)*dat$x$s.iid)/sqrt(tau.s),
                     dat$x$s.str)
    if (type==1) 
        ev.st <- list(values=rep(1, m*n), vectors=diag(m*n)) 
    if (type==2)
        ev.st <- list(values=rep(ev.t$values, each=n),
                      vectors=kronecker(ev.t$vectors, diag(n)))
    if (type==3) 
        ev.st <- list(values=rep(ev.s$values, m), 
                      vectors=kronecker(diag(m), ev.s$vectors))
    if (type==4) 
        ev.st <- list(values=rep(ev.t$values, each=n)*rep(ev.s$values, m),  
                      vectors=kronecker(ev.t$vectors, ev.s$vectors))
    dat$x$spacetime <- qsample(ev=ev.st)/sqrt(tau.st)
    dat$x$eta <- intercept + dat$x$time[dat$time] +
        dat$x$space[dat$space] + dat$x$spacetime[dat$spacetime]
    return(dat)
}

