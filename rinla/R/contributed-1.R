## Export: inla.contrib.sd

##!\name{contrib.sd}
##!\alias{contrib.sd}
##!\alias{inla.contrib.sd}
##!\title{Computes the standard deviation for the structured (random) effects in an INLA model}
##!\description{Computes the posterior distribution of the standard deviations for the structured 
##!             (random) effects in an INLA model, starting from the default output based on the precisions}
##!\usage{
##!    inla.contrib.sd(model, nsamples=1000)
##!}
##!\arguments{
##!    \item{model}{An INLA model, fitted calling the \code{inla()}-function. The formula specified for the model 
##!            should include at least one structured (random) effect in the form \code{f(variable, model="iid")}.}
##!    \item{nsamples}{The number of simulations from the posterior distribution of the standard deviations
##!                    used to compute the summary statistics}
##!}
##!\value{
##!     \code{inla.contrib.sd} returns a matrix \code{samples} including the simulated values from the
##!                            posterior distributions as well as a summary table \code{hyper} reporting 
##!                            the mean, sd and 95\% credible interval for the posterior distributions of each random effect.
##!}
##!\author{Gianluca Baio \email{gianluca@stats.ucl.ac.uk}}
##!\seealso{
##!    \code{\link{inla}}
##!}
##!\examples{
##!# Data generation
##!n=12
##!Ntrials = sample(c(80:100), size=n, replace=TRUE)
##!eta = rnorm(n,0,0.5)
##!prob = exp(eta)/(1 + exp(eta))
##!y = rbinom(n, size=Ntrials, prob = prob)
##!data=data.frame(y=y,z=1:n)
##!formula=y~f(z,model="iid")
##!m=inla(formula,data=data,family="binomial",Ntrials=Ntrials)
##!summary(m)
##!s=inla.contrib.sd(m)
##!s$hyper
##!hist(s$samples,xlab="standard deviation for z",main="")
##!}

`inla.contrib.sd` <- function(model, nsamples=1000)
{
    ## contributed by Gianluca Baio <gianluca@stats.ucl.ac.uk>

    ## Computes the sd for the random effects in an INLA model
    ## 1. Defines the precision (generates a matrix with bins and
    ## density on the precision scale)

    ## 2. Simulates replications from the posterior distributions of
    ## the quantities of interest

    ## Names of the variables associated with structured effects
    rand.effs <- names(model$marginals.hyperpar)
    for (i in 1:length(rand.effs)) {
        cmd <- paste("prec.marg.",i,"<-model$marginals.hyperpar$'",rand.effs[i],"'",sep="")
        eval(parse(text=cmd)) # marginal distribution of the precision, tau
        ## Simulation from the posterior marginal distribution for sigma = 1/sqrt(tau)
        cmd <- paste("sigma.", i,
                     "<- inla.rmarginal(nsamples,inla.marginal.transform(function(x) 1/sqrt(x), prec.marg.",
                     i,"))",sep="")
        eval(parse(text=cmd))
    }

    ## Outputs of the function
    mat <- matrix(NA, nsamples, length(rand.effs))
    for (i in 1:length(rand.effs)) {
        cmd <- paste("mat[,i] <- sigma.",i,sep="")
        eval(parse(text=cmd)) 
    }
    names2 <- gsub("Precision","sd",rand.effs)
    colnames(mat) <- names2

    tab <- matrix(NA,length(rand.effs),4)
    for (i in 1:length(rand.effs)) {
        tab[i,] <- c(mean(mat[,i]),sd(mat[,i]),quantile(mat[,i],.025),quantile(mat[,i],.975))
    }
    rownames(tab) <- names2
    colnames(tab) <- c("mean","sd","2.5%","97.5%")

    return (list(samples=mat, hyper=tab))
}
