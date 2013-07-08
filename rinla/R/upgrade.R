### The upgrade utility

`inla.update` = function(...)
{
    inla.upgrade(...)
}

`inla.upgrade` = function(lib = NULL, testing = FALSE, force = FALSE)
{
    ## include depends-on packages here
    for(p in c("sp", "Matrix", "splines")) {
        if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
            if (!require(p, quietly = TRUE, lib.loc = lib, character.only=TRUE)) {
                install.packages(p)
                ##stop(paste("INLA need package `", p, "' to be fully functional; please install", sep=""))
            }
        }
    }

    if (testing)
        www = "http://www.math.ntnu.no/inla/binaries/testing"
    else 
        www = "http://www.math.ntnu.no/inla/binaries"

    b.date = scan(paste(www,"/build.date", sep=""), quiet=TRUE, what = character(0))
    if (exists("inla.version")) {
        bb.date = inla.version("bdate")
    } else {
        bb.date = "INLA.is.not.installed"
    }

    if (b.date == as.character(bb.date)) {
        cat("\nYou have the newest version of INLA:\n")
        inla.version()
        if (!force)
            return (invisible())
        else
            cat("\nForce a new install\n")
    }
    
    ## download and install INLA
    if (inla.os("windows")) {
        suff = ".zip"
        tp = "win.binary"
    } else {    
        suff = ".tgz"
        tp = "source"
    }
    dfile = paste(tempdir(), .Platform$file.sep, "INLA", suff, sep="")
    sfile = paste(www, "/INLA", suff, sep="")
    download.file(sfile, dfile)


    ## use previous path if available
    if (is.null(lib)) {
        lib = searchpaths()[grep("[\\/]INLA$", searchpaths())]
        if (length(lib) == 0) {
            lib = NULL
        } else {
            while(length(grep("/?INLA$", lib)))
                lib = sub("/?INLA$", "", lib)
        }
        
        if (is.null(lib)) {
            ## ###########################################
            ## this part is copied from install.packages()
            ## ###########################################
            if (missing(lib) || is.null(lib)) {
                lib <- .libPaths()[1L]
                if (length(.libPaths()) > 1L) 
                    warning(gettextf("argument 'lib' is missing: using '%s'", 
                                     lib), immediate. = TRUE, domain = NA)
            }
            ok <- file.info(lib)$isdir & (file.access(lib, 2) == 0)
            if (length(lib) > 1 && any(!ok)) 
                stop(sprintf(ngettext(sum(!ok), "'lib' element '%s'  is not a writable directory", 
                                      "'lib' elements '%s' are not writable directories"), 
                             paste(lib[!ok], collapse = ", ")), domain = NA)
            if (length(lib) == 1 && inla.os("windows")) {
                ok <- file.info(lib)$isdir
                if (ok) {
                    fn <- file.path(lib, "_test_dir_")
                    unlink(fn, recursive = TRUE)
                    res <- try(dir.create(fn, showWarnings = FALSE))
                    if (inherits(res, "try-error") || !res) {
                        ok <- FALSE
                    } else {
                        unlink(fn, recursive = TRUE)
                    }
                }
            }
            if (length(lib) == 1L && !ok) {
                warning(gettextf("'lib = \"%s\"' is not writable", lib), 
                        domain = NA, immediate. = TRUE)
                userdir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"), 
                                           .Platform$path.sep))[1L]
                if (interactive() && !file.exists(userdir)) {
                    msg <- gettext("Would you like to create a personal library\n'%s'\nto install packages into?")
                    if (inla.os("windows")) {
                        ans <- winDialog("yesno", sprintf(msg, userdir))
                        if (ans != "YES") 
                            stop("unable to install the INLA package")
                    } else {
                        ans <- readline(paste(sprintf(msg, userdir), 
                                              " (y/n) "))
                        if (substr(ans, 1L, 1L) == "n") 
                            stop("unable to install the INLA package")
                    }
                    if (!dir.create(userdir, recursive = TRUE))
                        stop("unable to create ", sQuote(userdir))
                    lib <- userdir
                    .libPaths(c(userdir, .libPaths()))
                } else {
                    stop("unable to install packages")
                }
            }
            ## ###########################################
            ## end of copy...
            ## ###########################################
        }
    } else {
        ## ###########################################
        ## same here
        ## ###########################################
        ok <- file.info(lib)$isdir & (file.access(lib, 2) == 0)
        if (length(lib) > 1 && any(!ok)) 
            stop(sprintf(ngettext(sum(!ok), "'lib' element '%s'  is not a writable directory", 
                                  "'lib' elements '%s' are not writable directories"), 
                         paste(lib[!ok], collapse = ", ")), domain = NA)
    }

    ## remove old library before installing the new one
    try(detach(package:INLA), silent = TRUE)
    try(unloadNamespace("INLA"), silent = TRUE)

    install.packages(dfile, lib = lib, repos=NULL, type = tp)
    library(INLA, lib.loc = lib)

    cat("\nType\n\tinla.version()\nto display the new version of R-INLA. Thanks for upgrading.\n\n")
    return (invisible())
}
