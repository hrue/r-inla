## Export: inla.nmix.lambda.fitted


### The inla.nmix.lambda.fitted() function is a helper function for taking an 'nmix' or 'nmixnb'
### model, and computing expected abundance values for each site or site-by-year combination in
### the data, using the linear predictor for ## lambda (hereafter, 'fitted values').
###
### The uncertainty associated with fitted values derives from repeated sampling of INLA
### posteriors for the parameters of the linear predictor, and repeated solving of the linear
### predictor equation for each site or site-by-year combination. By default, fitted values are
### exponentiated to the count scale, and summaries of approximate fitted posteriors are returned
### by the function. Full estimated posteriors of fitted values are also available as an optional
### output.
###
### Tim Meehan <tmeehan@audubon.org>, 16 February 2018. 

##!\name{inla.nmix.lambda.fitted}
##!\alias{inla.nmix.lambda.fitted}
##!\alias{nmix.lambda.fitted}
##!\alias{inla.nmix.fitted}
##!
##!\title{
##!  Estimate posterior distributions of fitted lambda values
##!}
##!
##!\description{For use with \code{'nmix'} and \code{'nmixnb'} models. This function takes the
##!    information contained in an object returned by \code{inla()} and uses the contents to create
##!    fitted lambda values using the linear predictor for log(lambda), the input covariate values,
##!    and samples from the posteriors of the model hyperparameters. Fitted values from the linear
##!    predictor are exponentiated, by default, before being returned.}
##!
##!\usage{
##!inla.nmix.lambda.fitted(result, sample.size = 1000,
##!                        return.posteriors = FALSE, scale = "exp")
##!}
##!
##!\arguments{
##!    \item{result}{The output object from a call to \code{inla()}, where the family argument has
##!  been set to \code{'nmix'} or \code{'nmixnb'}. For the function to work, the call to
##!  \code{inla()} should also include the argument \code{control.compute=list(config = TRUE))}.}
##!
##!  \item{sample.size}{The size of the sample from the posteriors of the model hyperparameters.
##!  This sample size ends up being the size of the estimated posterior for a fitted lambda value.
##!  Default is 1000. Larger values are recommended.}
##!
##!  \item{return.posterior}{A logical value for whether or not to return the full estimated
##!  posteriors for each fitted value (\code{TRUE}), or just a summary of the posteriors
##!  (\code{FALSE}). Default is \code{FALSE}.}
##!
##!  \item{scale}{A character string, where the default string, \code{"exp"}, causes values from
##!  the linear predictor to be exponentiated before being returned. The string, \code{"log"},
##!  causes values to be returned on the \code{log(lambda)} scale.} 
##!}
##!
##!\value{
##!    \item{fitted.summary}{A data frame with summaries of estimated posteriors of fitted lambda
##!  values. The number of rows equals the number of rows in the data used to create the
##!  \code{'nmix'} or \code{'nmixnb'} model. There are six columns of summary statistics for each
##!  estimated posterior. Columns include an \code{index}, \code{mean.lambda}, \code{sd.lambda},
##!  \code{quant025.lambda}, \code{median.lambda}, \code{quant975.lambda}, and \code{mode.lambda}.}
##!
##!  \item{fitted.posteriors}{A data frame containing samples that comprise the full estimated
##!  posteriors of fitted values. The number of rows equals the number of rows in the data used to
##!  create the \code{'nmix'} or \code{'nmixnb'} model. The number of columns equals one plus the
##!  number of samples specified by the \code{sample.size} argument.} 
##!}
##!
##!\references{
##!    See documentation for families "nmix" and "nmixmb": \code{inla.doc("nmix")}
##!}
##!\author{
##!    Tim Meehan <tmeehan@audubon.org>
##!}
##!\note{
##!    This function is experimental.
##!}
##!
##!\examples{
##!## an example analysis of an N-mixture model using simulated data
##!## set parameters
##!n <- 75                       # number of study sites
##!nrep.max <- 5                 # number of surveys per site
##!b0 <- 0.5                     # lambda intercept, expected abundance
##!b1 <- 2.0                     # effect of x1 on lambda
##!a0 <- 1.0                     # p intercept, detection probability
##!a2 <- 0.5                     # effect of x2 on p
##!size <- 3.0                   # size of theta
##!overdispersion <- 1 / size    # for negative binomial distribution
##!
##!## make empty vectors and matrix
##!x1 <- c(); x2 <- c()
##!lambdas <- c(); Ns <- c()
##!y <- matrix(NA, n, nrep.max)
##!
##!## fill vectors and matrix
##!for(i in 1:n) {
##!    x1.i <- runif(1) - 0.5
##!    lambda <- exp(b0 + b1 * x1.i)
##!    N <- rnbinom(1, mu = lambda, size = size)
##!    x2.i <- runif(1) - 0.5
##!    eta <- a0 + a2 * x2.i
##!    p <- exp(eta) / (exp(eta) + 1)
##!    nr <- sample(1:nrep.max, 1)
##!    y[i, 1:nr] <- rbinom(nr, size = N, prob = p)
##!    x1 <- c(x1, x1.i); x2 <- c(x2, x2.i)
##!    lambdas <- c(lambdas, lambda); Ns <- c(Ns, N)
##!}
##!
##!## bundle counts, lambda intercept, and lambda covariates
##!Y <- inla.mdata(y, 1, x1)
##!
##!## run inla and summarize output
##!result <- inla(Y ~ 1 + x2,
##!  data = list(Y=Y, x2=x2),
##!  family = "nmixnb",
##!  control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01,
##!                      prec.intercept = 0.01),
##!  control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
##!                                    theta2 = list(param = c(0, 0.01)),
##!                                    theta3 = list(prior = "flat",
##!                                                 param = numeric()))),
##!  control.compute=list(config = TRUE)) # important argument
##!summary(result)
##!
##!## get and evaluate fitted values
##!lam.fits <- inla.nmix.lambda.fitted(result, 5000)$fitted.summary
##!plot(lam.fits$median.lambda, lambdas)
##!round(sum(lam.fits$median.lambda), 0); sum(Ns)
##!}


inla.nmix.lambda.fitted <- function(result, sample.size=1000,
                                    return.posteriors=FALSE,
                                    scale="exp")
{
    fam <- result$.args$family
    if (length(grep(pattern = "nmix", x = fam)) == 0) {
        stop("This function is only for models with 'nmix' or 'nmixnb' likelihoods")
    }
    if (missing(result)) 
        stop("Please specify a model result")
    s.check <- as.numeric(scale == "exp") + as.numeric(scale == "log")
    if (s.check == 0) 
        stop("Scale must be set to 'exp' or 'log'")
    if (sample.size < 500) 
        warning("Please increase the sample size")

    ## Get counts and lambda covariates from 'inla.mdata' object
    mdata.obj <- result$.args$data[[1]]
    counts <- as.data.frame(mdata.obj[grep(pattern="Y", names(mdata.obj))])
    lambda.covs <- as.data.frame(mdata.obj[grep(pattern="X", names(mdata.obj))])
    lambda.covs <- as.matrix(lambda.covs)
    n.counts <- ncol(counts)
    n.lambda.covs <- ncol(lambda.covs)
    n.data <- nrow(counts)

    ## Get samples from hyperpars
    hyperpar.samples <- inla.hyperpar.sample(sample.size, result)
    s.names <- rownames(hyperpar.samples)

    ## Discard overdispersion marginal posterior if 'nmixnb'
    if (fam == "nmixnb"){
        hyperpar.samples <- hyperpar.samples[,-(ncol(hyperpar.samples))]
    }
    n.samp.covs <- ncol(hyperpar.samples)
    if (n.lambda.covs != n.samp.covs) {
        stop("The number of hyperparameters and covariates does not match")
    }

    ## Combine lambda covariates and hyperparameter posteriors
    fitted.posteriors <- matrix(-1.0000, nrow=n.data, ncol=sample.size)
    ## For each site or site-by-year combination
    for(i in 1:n.data){
        obs <- lambda.covs[i,]
        ## For each sample from the hyperparameter posterior
        for(j in 1:sample.size){
            post <- hyperpar.samples[j, ]
            fitted <- sum(obs * post)
            fitted.posteriors[i,j] <- fitted
        }
    }

    ## Clean up the resulting matrix
    index <- 1:n.data
    fitted.posteriors <- as.data.frame(fitted.posteriors)
    row.names(fitted.posteriors) <- NULL
    names(fitted.posteriors) <- s.names

    ## Create posterior summaries of fitted values
    if (scale=="exp"){
        fitted.posteriors <- exp(fitted.posteriors)
    }
    fitted.meds <- round(apply(fitted.posteriors, 1, median), 4)
    fitted.means <- round(apply(fitted.posteriors, 1, mean), 4)
    fitted.sds <- round(apply(fitted.posteriors, 1, sd), 4)
    fitted.q025 <- round(apply(fitted.posteriors, 1, quantile, probs=0.025), 4)
    fitted.q500 <- round(apply(fitted.posteriors, 1, quantile, probs=0.500), 4)
    fitted.q975 <- round(apply(fitted.posteriors, 1, quantile, probs=0.975), 4)
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    fitted.modes <- round(apply(fitted.posteriors, 1, Mode), 4)
    fitted.summary <- data.frame(mean.lambda=fitted.means, sd.lambda=fitted.sds,
                                 quant025.lambda=fitted.q025,
                                 median.lambda=fitted.q500,
                                 quant975.lambda=fitted.q975,
                                 mode.lambda=fitted.modes)
    fitted.summary <- cbind(index, fitted.summary)
    fitted.posteriors <- cbind(index, fitted.posteriors)

    ## Create returned object
    if (return.posteriors==TRUE){
        out <- list(fitted.summary=fitted.summary,
                    fitted.posteriors=fitted.posteriors)
    } else {
        out <- list(fitted.summary=fitted.summary)
    }
    return(out)
}


