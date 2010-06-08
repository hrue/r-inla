## create a spde model using fmesher.

`inla.create.spde` = function(data.locations = list(x=NULL,y=NULL,z=NULL),
        boundary = list(x=NULL, y=NULL, z=NULL), dir = tempdir())
{

    ## create files using fmesher

    ## create the directory if not already there, otherwise, create a
    ## random one if dir=NULL.
    if (is.na(file.info(dir)$isdir) || !file.info(dir)$isdir) {
        if (!mkdirs(dir)) {
            stop(paste("Fail to create directory:", dir))
        }
    }

    if (!is.list(data.locations))
        stop("Argument `data.locations' is required.")

    n = length(data.locations$x)
    stopifnot(length(data.locations$y) == n)
    stopifnot(length(data.locations$z) == 0)

    ## this does not work for dimension 3.... FIXME!!! locations must
    ## be unique. check ?? make the locations complex numbers and then
    ## check...
    stopifnot(length(data.locations$x + 1i * data.locations$y) == n)

    ## Following FL's convension, the location file MUST end with a
    ## `s0'...  which means that we have to separate the location file
    ## and the argument to fmesher.
    A = cbind(data.locations$x, data.locations$y)
    loc.file.argument = paste(tempfile(tmpdir=dir), ".", sep="")
    loc.file = inla.write.fmesher.file(A, filename = paste(loc.file.argument, "s0", sep=""))

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
    
    return (loc.file.argument)
}
