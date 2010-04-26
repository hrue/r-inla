
`inla.os` = function(type = c("linux", "mac", "windows", "else"))
{
    if (missing(type))
        stop("Type of OS is required.")
    type = match.arg(type)
    
    if (type == "windows")
        return (.Platform$OS.type == "windows")
    else if (type == "mac")
        return (length(grep("mac", .Platform$pkgType)) > 0)
    else if (type == "linux")
        return ((.Platform$OS.type == "unix") && !inla.os("mac"))
    else if (type == "else")
        return (TRUE)
    else
        stop("This shouldn't happen.")
}
`inla.os.type` = function()
{
    for (os in c("windows", "mac", "linux", "else"))
        if (inla.os(os))
            return (os)
    stop("This shouldn't happen.")
}

