## Nothing to export

`inla.cygwin.check.path` = function(path = inla.getOption("cygwin"))
{
    return (file.exists(path) && file.info(path)$isdir)
}

`inla.cygwin.protect.special.chars` = function(cmd)
{
    ##cmd = gsub(" ", "\\\\ ", cmd)
    return (cmd) 
}

`inla.cygwin.run.command` = function(command, file.log = NULL, ...)
{
    if (inla.cygwin.check.path()) {
        cmd = paste(inla.getOption("cygwin"),
                     "/bin/bash.exe -c ",
                     shQuote(paste("export PATH=/bin:/usr/bin:$PATH;",
                                   inla.cygwin.protect.special.chars(command),
                                   inla.ifelse(is.null(file.log), "", paste(" > ", file.log)), 
                                   sep=""), type="cmd"), sep="")
        system(cmd, ...)
    } else {
        stop(paste("Cannot find the CYGWIN installation:", inla.getOption("cygwin")))
    }
}

`inla.cygwin.map.filename` = function(filename, windows2cygwin = TRUE)
{
    if (windows2cygwin) {
        fnm = inla.paste(inla.cygwin.run.command(paste("cygpath -u ",
            inla.cygwin.protect.special.chars(filename)),
            intern = TRUE, ignore.stderr = TRUE))
    } else {
        fnm = inla.paste(inla.cygwin.run.command(paste("cygpath -w -m -s ",
            inla.cygwin.protect.special.chars(filename)),
            intern = TRUE, ignore.stderr = TRUE))
    }
    return (fnm)
}


        
