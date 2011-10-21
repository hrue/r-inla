## Call fmesher

`inla.fmesher.call` = function(
        fmesher.call = inla.getOption("fmesher.call"),
        all.args, prefix)
{
    if (inla.os("linux") || inla.os("mac")) {
        echoc = system(paste(shQuote(fmesher.call), all.args, shQuote(prefix)))
    }
    else if(inla.os("windows")) {
        if (TRUE) {
            echoc = try(system(paste(shQuote(fmesher.call),
                                     all.args,
                                     shQuote(prefix))), silent=TRUE)
            echoc = 0
        } else {
            ## we might need it if we want one day to make the fmesher program run remotely
            stop("This code is not in use.")
            ##
            echoc = try(inla.cygwin.run.command(
                    paste(inla.cygwin.map.filename(fmesher.call),
                          all.args,
                          inla.cygwin.map.filename(prefix))), silent=TRUE)
            echoc = 0
        }
    }
    else
        stop("\n\tNot supported architecture.")
}

`inla.fmesher.make.prefix` = function(dir = tempdir(), prefix = NULL)
{
    if (is.null(prefix)) {
        if (is.null(dir))
            prefix = paste(tempfile(), ".", sep="")
        else {
            inla.dir.create(dir)
            prefix = paste(tempfile(tmpdir=dir), ".", sep="")
        }
    } else {
        inla.dir.create(inla.dirname(prefix))
    }
    return (prefix)
}




## Obsolete functions.

`inla.create.spde` = function(mesh=NULL, prefix = NULL,
        fem=2, sph0=NULL, sph=NULL, bspline=NULL, points2mesh=NULL)
{
    stop("'inla.create.spde' has been replaced by 'inla.spde.create'.")
}

`inla.fmesher.mesh` = function(
        locations,
        dir = tempdir(),
        prefix = NULL,
        rcdt = c(21,-1.0,-0.5),
        cet = c(8,-0.1),
        boundary = NULL,
        interior = NULL,
        x11.delay = NULL,
        min.input.distance = 0)
{
    stop("'inla.fmesher.mesh' has been replaced by 'inla.mesh.create'.")
}
