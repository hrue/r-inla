## create a spde model using fmesher.

`inla.create.spde` = function(locations, boundary = NULL, dir = tempdir())
{

    ## create files using fmesher

    ## create the directory if not already there, otherwise, create a
    ## random one if dir=NULL.
    if (is.na(file.info(dir)$isdir) || !file.info(dir)$isdir) {
        if (!dir.create(dir, recursive=TRUE)) {
            stop(paste("Fail to create directory:", dir))
        }
    }

    Loc = as.matrix(locations)
    n = dim(Loc)[1]
    ss = inla.unique.rows(Loc)
    s = ss$rows
    idx = ss$idx

    ## Following FL's convension, the location file MUST end with a
    ## `s0'...  which means that we have to separate the location file
    ## and the argument to fmesher.
    loc.file.argument = paste(tempfile(tmpdir=dir), ".", sep="")
    loc.file = inla.write.fmesher.file(s, filename = paste(loc.file.argument, "s0", sep=""))

    ## additional argumets
    all.args = paste("--rcdt", inla.getOption("fmesher.arg"))
    fmesher.call = inla.getOption("fmesher.call")
    
    if (inla.os("linux") || inla.os("mac")) {
        echoc = system(paste(shQuote(fmesher.call), all.args, loc.file.argument))
    }
    else if(inla.os("windows")) {
        if (!remote) {
            echoc = try(system(paste(shQuote(fmesher.call), all.args, loc.file.argument)), silent=TRUE)
            echoc = 0
        } else {
            echoc = try(inla.cygwin.run.command(
                    paste(inla.cygwin.map.filename(fmesher.call),
                          all.args,
                          inla.cygwin.map.filename(loc.file.argument))), silent=TRUE)
            echoc = 0
        }
    }
    else
        stop("\n\tNot supported architecture.")
    
    return (list(prefix = loc.file.argument, locations = s, locations.idx = idx))
}

