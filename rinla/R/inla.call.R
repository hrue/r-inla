`inla.call.builtin` <- function() {
    ## cannot call inla.getOption() here as it leads to an infinite recursive call. do this
    ## manually instead.
    opt.default <- inla.getOption.default()
    if (exists("inla.options", envir = inla.get.inlaEnv())) {
        opt <- get("inla.options", envir = inla.get.inlaEnv())
    } 

    if (inla.os("mac")) {
        fnm <- system.file(paste("bin/mac/", inla.os.32or64bit(), "bit/inla.mkl.run", sep = ""), package = "INLA")
    } else if (inla.os("mac.arm64")) {
        fnm <- system.file("bin/mac.arm64/inla.run", package = "INLA")
    } else if (inla.os("linux")) {
        fnm <- system.file(paste("bin/linux/", inla.os.32or64bit(), "bit/inla.mkl.run", sep = ""), package = "INLA")
    } else if (inla.os("windows")) {
        fnm <- system.file(paste("bin/windows/", inla.os.32or64bit(), "bit/inla.exe", sep = ""), package = "INLA")
    } else {
        stop("Unknown OS")
    }

    if (file.exists(fnm)) {
        return(fnm)
    } else {
        stop(paste("INLA installation error; no such file", fnm))
    }
}

`inla.call.no.remote` <- function() {
    ## return what is defined in options$inla.call except for 'remote', for which we revert back
    ## to the builtin one
    inla.call <- inla.getOption("inla.call")
    if (is.null(inla.call) || any(inla.strcasecmp(inla.call, c("remote", "inla.remote")))) {
        inla.call <- inla.call.builtin()
    }
    return(inla.call)
}

`inla.fmesher.call.builtin` <- function() {
    if (inla.os("mac")) {
        fnm <- system.file(paste("bin/mac/", inla.os.32or64bit(), "bit/fmesher.run", sep = ""),
                           package = "INLA")
    } else if (inla.os("mac.arm64")) {
        fnm <- system.file("bin/mac.arm64/fmesher.run", package = "INLA")
    } else if (inla.os("linux")) {
        fnm <- system.file(paste("bin/linux/", inla.os.32or64bit(), "bit/fmesher.run", sep = ""), package = "INLA")
    } else if (inla.os("windows")) {
        fnm <- system.file(paste("bin/windows/", inla.os.32or64bit(), "bit/fmesher.exe", sep = ""), package = "INLA")
    } else {
        stop("Unknown OS")
    }

    if (file.exists(fnm)) {
        return(fnm)
    } else {
        stop(paste("INLA installation error; no such file", fnm))
    }
}
