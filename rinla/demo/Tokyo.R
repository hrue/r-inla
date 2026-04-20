## Load the data
data(Tokyo)
summary(Tokyo)

## Define the model
formula = y ~ f(time, model="rw2", cyclic=TRUE, param=c(1,0.0001)) - 1

## The call to inla
result = inla(formula, family="binomial", Ntrials=n, data=Tokyo)

## Once more: the call to inla in verbose mode
result = inla(formula, family="binomial", Ntrials=n, data=Tokyo, verbose = TRUE)

## Summarise the results
summary(result)

## Plot the results
plot(result)

if (FALSE) {
    ## only for classic mode
    h = inla.hyperpar(result)
    summary(h)
    plot(h)
}

