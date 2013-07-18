## Export:  inla.ssh.copy.id

##!\name{inla.ssh.copy.id}
##!\alias{inla.ssh.copy.id}
##!\alias{ssh.copy.id}
##!
##!\title{Setup remote-access}
##!
##!\description{Print the path to the internal script to transfer ssh-keys}
##!
##!\usage{
##!inla.ssh.copy.id()
##!}
##!\arguments{
##!}
##!
##!\value{%%
##!  The path to the internal script to transfer ssh-keys
##!}
##!%%
##!
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!
##!\examples{
##!##TODO
##!}


`inla.ssh.copy.id` = function ()
{
    ## print the path to the copy-id script

    f = system.file("bin/remote/ssh-copy-id", package="INLA")
    if (inla.os("windows"))
        return (inla.cygwin.map.filename(f))
    else
        return (f)
}

`inla.remote` = function()
{
    if (inla.os("windows") && !inla.cygwin.check.path()) {
        cat("\n\n\tRemote computing in INLA from Windows require CYGWIN.\n")
        cat("\tThe CYGWIN-installation is not found. Please install it from www.cygwin.com\n")
        cat("\tMake sure to install also all `ssh' components; search for it, select them and install.\n")
        cat("\tIts easiest to use the standard installation path; c:/cygwin, otherwise give the new path using\n")
        cat("\t\tinla.setOption(\"cygwin\", \"NEWPATH\")\n")
        cat("\n")
        return (NULL)
    }

    inla=system.file("bin/remote/inla.remote", package="INLA")
    inla.setOption("inla.call", inla)

    f = paste(inla.get.HOME(), "/.inlarc", sep="")
    if (!file.exists(f)) {
        if (inla.os("windows"))
            ini=system.file("bin/remote/dot.inlarc-example-windows", package="INLA")
        else
            ini=system.file("bin/remote/dot.inlarc-example", package="INLA")

        if (!file.exists(ini))
            stop(paste("Missing file in installation:", ini, "This should not happen..."))
        file.copy(ini, f)
        cat(paste("Remote computing: Please edit remote-host details in\n\t\t", f, "\n"))
    } else {
        cat(paste("Remote computing is defined in\n\t\t", f, "\n"))
    }
    if (inla.os("windows")) {
        cat("\n\tIf you need to copy ssh-keys, you can do that using this script from within a cygwin-terminal\n")
        cat(paste("\t\t", inla.cygwin.map.filename(system.file("bin/remote/ssh-copy-id", package="INLA")), "\n\n"))
    }

    return (invisible())
}
