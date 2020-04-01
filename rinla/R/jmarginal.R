## Export: inla.rjmarginal 

##! \name{joint.marginal}
##! \alias{inla.joint.marginal}
##! \alias{joint.marginal}
##! \alias{rjmarginal}
##! \alias{inla.rjmarginal}
##!
##! \title{Sample from a joint marginal approximation}
##! 
##! \description{Sample from from a joint marginal approximation
##! as returned using argument \code{selection} in \code{inla}.
##! }
##! 
##! \usage{
##! inla.rjmarginal(n, jmarginal)
##! }
##! \arguments{
##!   \item{n}{The number of samples}
##!   \item{jmarginal}{A marginal object given  either by a \code{inla} object
##!     or \code{result$selection}}
##! }
##! 
##! \value{%%
##! A list with the samples in \code{samples} (matrix) and the corresponding log-densities
##! in \code{log.density} (vector).
##! Each column in \code{samples} contains one sample.
##! }
##! 
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \seealso{\code{\link{inla}}}
##!
##! \examples{
##! n = 10
##! x = 1+rnorm(n)
##! xx = 3 + rnorm(n)
##! y = 1 + x + xx + rnorm(n)
##! selection = list(xx=1, Predictor = 3:4, x=1)
##! r = inla(y ~ 1 + x + xx,
##!          data = data.frame(y, x, xx),
##!          selection = selection)
##! ns = 100
##! xx = inla.rjmarginal(ns, r)
##! 
##! print(cbind(mean = r$selection$mean, sample.mean = rowMeans(xx$samples)))
##! print("cov matrix")
##! print(round(r$selection$cov.matrix, dig=3))
##! print("sample cov matrix")
##! print(round(cov(t(xx$samples)), dig=3))
##! 
##! skew = function(z) mean((z-mean(z))^3)/var(z)^1.5
##! print(round(cbind(skew = r$selection$skewness,
##!                   sample.skew = apply(xx$sample, 1, skew)), dig=3))
##!}
##! 

`inla.rjmarginal` = function(n, jmarginal) 
{
    if (missing(jmarginal) || missing(n) || n <= 0) {
        return (list(samples = matrix(ncol=0, nrow=0), log.density = numeric(0)))
    }
    
    if (inherits(jmarginal, "inla")) {
        jmarginal = jmarginal$selection
    } else if (inherits(jmarginal, "inla.selection")) {
        ## ok
    } else {
        stop("Unknown object: argument 'jmarginal'")
    }

    ## build a fake 'inla'-object, so we can feed that into 'inla.posterior.sample'
    m = length(jmarginal$mean)
    r = list(
        misc = list(
            configs = list(
                n = m, 
                nconfig = 1,
                max.log.posterior = 0,
                nz = m*(m+1)/2, 
                contents = list(
                    tag = jmarginal$names,
                    start = 1:m,
                    length = rep(1, m)),
                ntheta = 0,
                config = list(list(theta = NULL,
                                   mean = jmarginal$mean,
                                   Q = inla.as.sparse(solve(jmarginal$cov.matrix)),
                                   log.posterior = 0,
                                   improved.mean = jmarginal$mean,
                                   Qinv = inla.as.sparse(jmarginal$cov.matrix),
                                   log.posterior.orig = 0,
                                   skewness = jmarginal$skewness)),
                constr = NULL
            )
        )
    )
    class(r) = "inla"

    x = inla.posterior.sample(n, r, use.improved.mean = TRUE, skew.corr = TRUE, add.names = FALSE)
    xx = matrix(unlist(lapply(x, function(z) z$latent)), ncol = n)
    log.dens = unlist(lapply(x, function(z) z$logdens$joint))
    rownames(xx) = jmarginal$names
    colnames(xx) = paste0("sample", 1:n)
    names(log.dens) = colnames(xx)
    
    return (list(samples = xx, log.density = log.dens))
}
