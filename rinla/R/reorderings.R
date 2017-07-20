## Export: inla.reorderings

##!\name{inla.reorderings}
##!\alias{inla.reorderings}
##!\alias{reorderings}
##!
##!\title{Reorderings methods for sparse matrices}
##!
##!\description{Provide the names of all implemented reordering schemes}
##!
##!\usage{
##! inla.reorderings()
##!}
##!
##!\arguments{
##! None
##!}
##!
##!\value{The names of all available reorderings}
##!%%
##!
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!
##!\examples{
##!inla.reorderings()
##!}



## provide list of all reordering and map their names and code. OOPS:
## these mappings must be in sync with the ones in
## gmrflib/sparse-interface.h

`inla.reorderings.list` = function()
{
    r = list(
            auto = -1L,
            default = 0L,
            identity = 1L,
            band = 2L,
            metis = 3L,
            genmmd = 4L,
            amd = 5L,
            amdbar = 6L,
            md = 7L,
            mmd = 8L,
            amdc = 9L,
            amdbarc = 10L, 
            reverseidentity = 11L
            )
    return (r)
}

`inla.reorderings` = function()
{
    return (names(inla.reorderings.list()))
}

`inla.reorderings.name2code` = function(name)
{
    r = inla.reorderings.list()
    name = tolower(name)
    code = c()
    for(nam in name) {
        idx = which(nam == names(r))
        if (length(idx) > 0L) {
            code = c(code, as.integer(r[idx]))
        } else {
            stop(paste("Unknown reordering:", nam))
        }
    }
    return (code)
}

`inla.reorderings.code2name` = function(code)
{
    r = inla.reorderings.list()
    name = c()
    for(co in code) {
        idx = which(co == r)
        if (length(idx) > 0L) {
            name = c(name, names(r)[idx])
        } else {
            stop(paste("Unknown code:", co))
        }
    }
    return (name)
}
