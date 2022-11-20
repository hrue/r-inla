## Export: inla.likelihood

## ! \name{inla.likelihood}
## ! \alias{likelilhood}
## !
## ! \title{Compute log-likelihood and CDF for new data using the likeilhood from a inla-object}
## !
## ! \description{This function computes the log-likeilhood and computes the CDF for new data
## ! using the likelihood from a inla-object}
## ! \usage{
## !     inla.likelihood(y, type = c("loglik", "CDF"), linear.predictor = numeric(0), 
## !                     family = "gaussian", theta = NULL, E = 1, scale = 1,  
## !                     Ntrials = 1,  strata = 1,
## !                     link.model = "default", link.covariates = NULL)
## ! }
## !
## ! \arguments{
## !   \item{y}{Data point/vector}
## !   \item{type}{Return log-liklihood or the CDF?}
## !   \item{linear.predictor}{A vector of linear predictors for the liklihood/CDF to be
## !         evaluated at}
## !   \item{family}{The family}
## !   \item{theta}{Vector of hyperparameters in the internal scale, that goes into this likelihood}
## !   \item{E}{Constrant E for Poisson-like likehoods}
## !   \item{scale}{scaling, like for Gaussian likehood}
## !   \item{Ntrails}{Number of trials in Binomial-like likelihood}
## !   \item{link.model}{The link-model used}
## !   \item{link.covariates}{Possible covariates that goes into the link.model}
## !}
## !\details{Details goes here}
## !\value{value goes here}
## !\author{Havard Rue \email{hrue@r-inla.org}}
## !
## !\examples{
## ! }
 

`poisson.likelihood.d` <- function(y, linear.predictor, E, inv.link.function, ...)
{
    return (dpois(y, lambda = E * inv.link.function(linear.predictor), log = TRUE))
}

`poisson.likelihood.p` <- function(y, linear.predictor, E, inv.link.function, ...)
{
    return (ppois(y, lambda = E * inv.link.function(linear.predictor)))
}

`gaussian.likelihood.d` <- function(y, linear.predictor, scale, inv.link.function,
                                    theta, ...) {
    s <- 1 / sqrt(exp(theta) * scale)
    return (dnorm(y, m = inv.link.function(linear.predictor), sd = s, log = TRUE))
}

`gaussian.likelihood.p` <- function(y, linear.predictor, scale, inv.link.function,
                                    theta, ...) {
    s <- 1 / sqrt(exp(theta[1]) * scale)
    return (pnorm(y, m = inv.link.function(linear.predictor), sd = s))
}

`inla.surv.likelihood.d` <- function(...) {
    inla.surv.likelihood.core(..., internal.type = "d")
}

`inla.surv.likelihood.p` <- function(...) {
    inla.surv.likelihood.core(..., internal.type = "p")
}


## for weibull,  scale=exp(lin.predictor),  usually
weibull.likelihood.d <- function(y, alpha, scale) {
    return (dweibull(y, shape = alpha, scale = scale, log = TRUE))
}
weibull.likelihood.p <- function(y, alpha, scale) {
    return (pweibull(y, shape = alpha, scale = scale))
}

`inla.surv.likelihood.core` <- function(...) {
    args <- list(...)
    stopifnot(args$internal.type == "d")

    ## its possible to write this pretty generic later...
    if (args$family.arg.str$family == "weibull") {
        stopifnot(args$family.arg.str$variant == 1)

        ## should we use 'alpha' or 'theta' ? 

        y <- args$y.surv$time
        event <- args$y.surv$event
        shape <- args$family.arg.str$alpha
        scale <- args$inv.link.function(args$linear.predictor)
        truncation <- args$y.surv$truncation
        lower <- args$y.surv$lower
        upper <- args$y.surv$upper

        if (truncation > 0) {
            F.trunc <- weibull.likelihood.p(truncation, alpha = shape, scale = scale)
            FF.trunc <- 1.0/(1.0 - F.trunc)
        } else {
            F.trunc <- 0
            FF.trunc <- 1.0
        }

        if (event == 1 || event == 4) {
            ld <- weibull.likelihood.d(y, alpha = shape, scale = scale)
        } else {
            ld <- NA
        }

        if (event == 0 || event == 3 || event == 4) {
            F.lower <- weibull.likelihood.p(lower, alpha = shape, scale = scale)
            F.lower <- (F.lower - F.trunc) * FF.trunc
        } else {
            F.lower <- NA
        }

        if (event == 2 || event == 3 || event == 4) {
            F.upper <- weibull.likelihood.p(upper, alpha = shape, scale = scale)
            F.upper <-  (F.upper - F.trunc) * FF.trunc
        } else {
            F.upper <- NA
        }

        ## do we need to rebuild 'cure.prob' ? 
        pcure <- args$cure.prob

        ## this copied from the function loglikelihood_generic_surv_NEW in inla.c
        if (event == 1) {
            ## EVENT_FAILURE
            return (ld + log(FF.trunc))
        } else if (event == 0) {
            ## EVENT_RIGHT
            return (log(pcure + (1.0 - pcure) * (1.0 - F.lower)))
        } else if (event == 2) {
            ## EVENT_LEFT
            return (log((1.0 - pcure) * F.upper))
        } else if (event == 3) {
            ## EVENT_INTERVAL
            return (log((1.0 - pcure) * (F.upper - F.lower)))
        } else if (event == 4) {
            ## EVENT_ININTERVAL
            return (log(1.0 - pcure) + ld - log(F.upper - F.lower))
        } else {
            stop("This should not happen")
        }
    } else {
            stop("NOT YET IMPLEMENTED")
    }
}

`inla.likelihood` <- function(y = NULL, y.surv = NULL, type = c("loglik", "CDF"), linear.predictor = NULL, 
                              family = "gaussian", theta = NULL, E = 1, scale = 1,  
                              Ntrials = 1,  strata = 1,
                              cure.prob = 0, cure.beta = c(), cure.covariates = c(),
                              family.arg.str = NULL, 
                              link.model = NULL, link.covariates = NULL) 
{
    type <- match.arg(type)
    if (family == "normal") family <- "gaussian"
    
    ## for surv-models,  then the link is defined in family.arg.str$...
    if (!is.null(link.model)) {
        inv.link.function <- eval(parse(text = paste0("inla.link.inv", tolower(link.model))))
    } else {
        inv.link.function <- eval(parse(text = paste0("inla.link.inv",
                                                      tolower(family.arg.str$link.model))))
    }
    fun <- eval(parse(text = paste0(family, ".likelihood.", if (type == "loglik") "d" else "p")))

    return (fun(y = y,
                y.surv = y.surv,
                theta = theta,
                E = E,
                scale = scale, 
                Ntrials = Ntrials,
                strata = strata,
                inv.link.function = inv.link.function,
                link.covariates = link.covariates, linear.predictor = linear.predictor,
                ## survival stuff
                cure.prob = cure.prob,
                cure.beta = cure.beta,
                cure.covariates = cure.covariates,
                family.arg.str = family.arg.str))
}
