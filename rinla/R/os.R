## nothing to export


`inla.os` = function(type = c("linux", "mac", "windows", "else"))
{
    if (missing(type)) {
        stop("Type of OS is required.")
    }
    type = match.arg(type)
    
    if (type == "windows") {
        return (.Platform$OS.type == "windows")
    } else if (type == "mac") {
        result = (file.info("/Library")$isdir && file.info("/Applications")$isdir)
        if (is.na(result)) {
            result = FALSE
        }
        if (result) {
            ## check that the version is at least the one use to build the binaries. 
            s = system("sw_vers -productVersion", intern=T)
            vers = as.integer(strsplit(s, ".", fixed=TRUE)[[1]])
            ver = vers[1] + vers[2]/10
            s.req = 10.10 ## @@@HARDCODED@@@
            if (ver < s.req) {
                stop("Your version, ", s, ", of MacOSX is to old for R-INLA. Update MacOSX to at least version ",
                     as.character(s.req), sep="")
            }
        }
        return (result)
    } else if (type == "linux") {
        return ((.Platform$OS.type == "unix") && !inla.os("mac"))
    } else if (type == "else") {
        return (TRUE)
    } else {
        stop("This shouldn't happen.")
    }
}
`inla.os.type` = function()
{
    for (os in c("windows", "mac", "linux", "else")) {
        if (inla.os(os)) {
            return (os)
        }
    }
    stop("This shouldn't happen.")
}

## this seems be to a nice way to test 32/64 bits architecture.
`inla.os.32or64bit` = function()
{
    return (ifelse(.Machine$sizeof.pointer == 4, "32", "64"))
}
`inla.os.is.32bit` = function()
{
    return (inla.os.32or64bit() == "32")
}
`inla.os.is.64bit` = function()
{
    return (inla.os.32or64bit() == "64")
}
