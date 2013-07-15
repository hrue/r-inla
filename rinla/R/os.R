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
