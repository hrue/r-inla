## Nothing to Export.

`inla.dev.new` = function(...)
{
    ## If running in RStudio then don't open a new device,  otherwise,  do.
    dev = getOption("device")
    if (is.character(dev) && inla.strncasecmp(dev, "RStudioGD")) {
        return (NULL)
    } else {
        return (dev.new(...))
    }
}
