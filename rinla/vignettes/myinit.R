## Sourced by most vignettes as: if (file.exists("myinit.R")) source("myinit.R")
##
## Anything set here overrides what the session already worked out, so an
## unconditional inla.call sends every vignette to one developer's machine and
## the whole vignette build fails everywhere else with "No inla binary is
## installed". Keep the override, but only take it when that binary is actually
## there, so it applies on the machine it was written for and is inert on CI
## and on anyone else's checkout.
local({
    dev.inla <- "/home/hrue/bin/inla.mkl.work"
    if (file.exists(dev.inla)) {
        inla.setOption(inla.call = dev.inla)
        inla.setOption(num.threads = "8:1:2")
    } else {
        inla.setOption(num.threads = "2:1:1")
    }
})
