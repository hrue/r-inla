## RCSId = $Id: Tokyo-compare.R,v 1.7 2009/09/28 08:36:56 hrue Exp $
##
## this example show how to compare inla-results with MCMC using the
## buildt-in sampler using the simple Tokyo-example; see demo(Tokyo)

data(Tokyo)
formula = y ~ f(time, model="rw2", cyclic=TRUE, param=c(1,0.0001)) - 1

## we need to keep the results in a known place
wd.inla = tempfile()

## get the results using some important arguments..
result = inla(formula, family="binomial", Ntrials=n, data=Tokyo,
    ##
    ## where to place the results
    ##
    working.directory = wd.inla,
    ##
    ## we need to keep the results
    ##
    keep=TRUE,
    ##
    ## this is to prevent the results to be stored in binary-format
    ##
    inla.arg = "")

## same with mcmc-results...
wd.mcmc = tempfile()

## same as before
result = inla(formula, family="binomial", Ntrials=n, data=Tokyo,
    ##
    ## where to place the results
    ##
    working.directory = wd.mcmc,
    ##
    ## we need to keep the results
    ##
    keep = TRUE,
    ##
    ## secret options to inla: 5000 iterations using thinning of 10
    ## and the scaling of the hyperparmeter proposal is 1.
    ##
    inla.arg = "-m mcmc -N 5000 -T 10 -S 1")

## enter the results to the comparison-routine, which interactively
## let you compare the results
inla.compare.results(wd.inla, wd.mcmc)

## remove the working-files
unlink(wd.inla, recursive = TRUE)
unlink(wd.mcmc, recursive = TRUE)
