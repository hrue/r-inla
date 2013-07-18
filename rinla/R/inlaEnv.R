## Nothing to export

## Define the environment for INLA used to store options and the
## model-list. Reuse the environment if it is there already.

if (exists(".inlaEnv") && is.environment(.inlaEnv)) {
    ## then reuse it
} else {
    .inlaEnv = new.env()
}

`inla.get.inlaEnv` = function(...)
{
    if (exists(".inlaEnv") && is.environment(.inlaEnv))
        return (.inlaEnv)
    stop("Environment '.inlaEnv' does not exists and is required for INLA to work. Restart 'R'.")
}

