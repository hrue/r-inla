## Export: summary!inla.jmarginal
## Export: print!summary.inla.jmarginal
## Export: print!inla.jmarginal

## Export: inla.rjmarginal inla.rjmarginal.eval

##! \name{joint.marginal}
##! \alias{inla.joint.marginal}
##! \alias{inla.joint.marginal.eval}
##! \alias{joint.marginal}
##! \alias{joint.marginal.eval}
##! \alias{rjmarginal}
##! \alias{rjmarginal.eval}
##! \alias{inla.rjmarginal}
##! \alias{inla.rjmarginal.eval}
##!
##! \title{Sample and evaluate from a joint marginal approximation}
##! 
##! \description{Sample and evalue from from a joint marginal approximation
##! as returned using argument \code{selection} in \code{inla}.
##! }
##! 
##! \usage{
##! inla.rjmarginal(n, jmarginal)
##! inla.rjmarginal.eval(fun, samples, ...) 
##! }
##! \arguments{
##!   \item{n}{The number of samples}
##!   \item{jmarginal}{A marginal object given  either by a \code{inla} object
##!     or \code{result$selection}}
##!   \item{fun}{A function which is evaluated for each sample, similar to
##!              \code{inla.posterior.sample.eval}: please see the documentation
##!              for this functions for details. }
##!   \item{samples}{The samples, as in the form of the output from \code{inla.rjmarginal}}
##! }
##! 
##! \value{%%
##! \code{inla.rjmarginal} returns a list with the samples in \code{samples}
##! (matrix) and the corresponding log-densities
##! in \code{log.density} (vector). Each column in \code{samples} contains one sample.
##!
##! \code{inla.rjmarginal.eval} returns a matrix, where each row is the (vector) function
##! evaluated at each sample.
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
##!
##! ## illustrating the eval function
##! n = 10
##! x = rnorm(n)
##! eta = 1 + x
##! y = eta + rnorm(n, sd=0.1)
##! selection = list(x = 1, Predictor = c(1, 2, 4, 5),  '(Intercept)' = 1)
##! r = inla(y ~ 1 + x,
##!          data = data.frame(y, x),
##!          selection = selection)
##! xx = inla.rjmarginal(100,  r)
##! xx.eval = inla.rjmarginal.eval(function() c(x, Predictor, Intercept),  xx)
##! print(cbind(xx$samples[, 1]))
##! print(cbind(xx.eval[, 1]))
##!}


`inla.rjmarginal` = function(n, jmarginal) 
{
    if (missing(jmarginal) || missing(n) || n <= 0) {
        return (list(samples = matrix(ncol=0, nrow=0), log.density = numeric(0)))
    }
    
    if (inherits(jmarginal, "inla")) {
        jmarginal = jmarginal$selection
    } else if (inherits(jmarginal, "inla.jmarginal")) {
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
    colnames(xx) = paste0("sample:", 1:n)
    names(log.dens) = colnames(xx)
    
    return (list(samples = xx, log.density = log.dens))
}

`inla.rjmarginal.eval` = function (fun, samples, ...) 
{
    stopifnot(all(names(samples) == c("samples", "log.density")))

    var = "inla.jmarginal.eval.warning.given"
    if (!(exists(var, envir = inla.get.inlaEnv()) &&
          get(var, envir = inla.get.inlaEnv()) == TRUE)) {
        warning("Function 'inla.rjmarginal.eval()' is experimental.")
        assign(var, TRUE, envir = inla.get.inlaEnv())
    }

    nm.split = strsplit(rownames(samples$samples), ":")
    ## determine the names
    nms = unique(unlist(lapply(nm.split, function(x) x[1])))
    contents = rep(list(list()), length(nms))
    for(i in seq_along(nms)) {
        contents[[i]]$tag = nms[i]
        contents[[i]]$start = min(which(nms[i] == lapply(nm.split, function(x) x[1])))
        contents[[i]]$sample.idx = which(nms[i] == lapply(nm.split, function(x) x[1]))
    }

    my.fun = function(a.sample, .contents, .fun, ...) {
        env = new.env()
        for (i in seq_along(.contents)) {
            .tmp = a.sample[contents[[i]]$sample.idx]
            assign(.contents[[i]]$tag, .tmp, envir = env)
        }
        if (exists("(Intercept)", envir = env)) {
            assign("Intercept", get("(Intercept)", envir = env), 
                   envir = env)
        }
        parent.env(env) = .GlobalEnv
        environment(.fun) = env
        return(.fun(...))
    }

    ret = apply(samples$samples, 2, my.fun, .fun = fun, .contents = contents, ...)
    ret = matrix(ret, ncol = ncol(samples$samples))
    colnames(ret) = paste0("sample:", 1:ncol(ret))
    rownames(ret) = paste0("fun[", 1:nrow(ret), "]")

    return(ret)
}

## simple methods to make some sensible outpout for r$selection, which also hide the .private
## section

`summary.inla.jmarginal` <- function(object) 
{
    ret <- c("Joint marginal is computed for:\n",
             paste("\t", object$names, "\n"))
    class(ret) <- "summary.inla.jmarginal"
    return (ret)
}

`print.summary.inla.jmarginal` <- function(x) 
{
    cat(x)
}

`print.inla.jmarginal`  <- function(x) 
{
    x$.private <- NULL
    class(x) <- class(list())
    print(x)
}
