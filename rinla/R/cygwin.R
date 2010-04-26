

`inla.cygwin.check.path` = function(path = inla.getOption("cygwin"))
{
    return (file.exists(path) && file.info(path)$isdir)
}

`inla.cygwin.run.command` = function(command, ...)
{
    if (inla.cygwin.check.path()) {
        cmd = paste(inla.getOption("cygwin"),
                     "/bin/bash.exe -c ",
                     shQuote(paste("export PATH=/bin:/usr/bin:$PATH;",
                                   command, sep=""), type="cmd"), sep="")
        system(cmd, ...)
    } else {
        stop(paste("Cannot find the CYGWIN installation:", inla.getOption("cygwin")))
    }
}

`inla.cygwin.map.filename` = function(filename, windows2cygwin = TRUE)
{
    if (windows2cygwin)
        return (inla.paste(inla.cygwin.run.command(paste("cygpath -u ", filename), intern = TRUE, ignore.stderr = TRUE)))
    else
        return (inla.paste(inla.cygwin.run.command(paste("cygpath -w -m -s ", filename), intern = TRUE, ignore.stderr = TRUE)))
}


        
