## copied from the mgcv-code...

## this does not work...
## invisible(utils::globalVariables(c("low", "high", "spde", "internal")))
## invisible(utils::suppressForeignCheck(c("low", "high", "spde", "internal")))

inla.print.version <- function() {
    info <- library(help = INLA)$info[[1]]
    if (!is.null(info)) {
        version <- info[pmatch("Version", info)]
        built <- strsplit(info[pmatch("Built", info)], "; ")[[1]][3]
        date <- info[pmatch("Date", info)]
        um <- strsplit(version, " ")[[1]]
        version <- um[nchar(um) > 0][2]
        um <- strsplit(date, " ")[[1]]
        date <- um[nchar(um) > 0]

        ## Label-and-value in a dotted column, the same layout inla.version()
        ## already uses, so the two agree. The previous form was four sentences
        ## in three different moods ("See ...", "List ... with ...",
        ## "Use ... to ..."), which has to be read rather than scanned, and it
        ## buried the two things a user actually types mid-sentence. Commands
        ## now sit at the end of the line, where the eye stops.
        ## Seconds are dropped from the build stamp: the date is what anyone
        ## checks, the clock time is noise in a greeting.
        built <- strsplit(built, " ")[[1]][1]
        pad <- function(x) {
            n <- 25L - nchar(x) - 1L
            paste0(x, " ", paste(rep(".", max(n, 1L)), collapse = ""))
        }
        hello <- paste0(
            "This is INLA_", version, ", built ", built, ".\n",
            "        ", pad("Help"),          ": www.r-inla.org/contact-us\n",
            "        ", pad("Models"),        ": inla.list.models()\n",
            "        ", pad("Documentation"), ": inla.doc(<NAME>)\n",
            "        ", pad("Powered by"),    ": sTiles, www.r-inla.org/sTiles"
        )
        ## Upgrade notice, disabled. It fetched
        ## https://inla.r-inla-download.org/VERSIONS on every attach and compared
        ## the running version against the published stable/testing ones. That
        ## server carries the upstream builds, so a version built here is never
        ## one of them and the notice fired every time, telling the user to
        ## "upgrade" to something older than what they are running.
        ## Wrapped rather than deleted so it can be restored in one edit.
        if (FALSE) {
            opts <- options()
            options(timeout = 2)
            suppressWarnings({
                vers <- try(readLines("https://inla.r-inla-download.org/VERSIONS",
                                      n = 4, encoding = "UTF-8"), silent = TRUE)
            })
            if (!inherits(vers, "try-error") && length(vers) == 4) {
                rem.space <- function(x) gsub("[ ]+","", x)
                stable <- rem.space(vers[1])
                testing <- rem.space(vers[2])
                major <- rem.space(vers[3])
                minor <- rem.space(vers[4])
                minor <- rem.space(strsplit(minor, "[.]")[[1]][1])
                current <- getNamespaceVersion("INLA")
                if (!(current == stable || current == testing)) {
                    majo <- rem.space(R.Version()$major)
                    mino <- rem.space(strsplit(R.Version()$minor, "[.]")[[1]][1])
                    if (majo != major || mino != minor) {
                        rstr <- paste0(" (require R-", major, ".", minor, ")")
                    } else {
                        rstr <- "."
                    }
                    hello <- paste0(hello, "\n",
                                    paste0(" - Consider upgrading R-INLA to testing[",  testing,
                                           "] or stable[", stable, "]", rstr))
                }
            }
            options(opts)
        }
        packageStartupMessage(hello)
    }
}

.onLoad <- function(...) {
    ## nothing for the moment
}
## No binary yet? Say so, once per session, and stop there.
##
## The package ships launcher scripts (inst/bin/<platform>/64bit/inla.run is
## ~2 KB, inla.mkl.run is 8 bytes), not a solver: that is a separate ~100 MB
## download. So a fresh install cannot run a model, and the first sign of it
## used to be an obscure failure inside inla(). This prints the one command
## that fixes it.
##
## It NEVER prompts and NEVER downloads on its own. library(INLA) must behave
## identically in a script, a container and a terminal, so nothing here can
## block on stdin or start a large transfer the caller did not ask for.
## Set options(inla.stiles.autoinstall = TRUE), or the environment variable
## INLA_STILES_AUTOINSTALL, to have it install unattended instead: that is
## opt-in precisely because it is the surprising behaviour.
##
## Existence is not the test. Those shipped launchers exist on every install,
## so file.exists() is always TRUE and would silence this permanently. Ask the
## binary whether it answers -V instead.
`inla.first.run.binary` <- function() {
    call <- tryCatch(inla.getOption("inla.call"), error = function(e) NULL)
    if (!is.null(call) && is.character(call) && nzchar(call) && file.exists(call)) {
        out <- suppressWarnings(tryCatch(
            system2(call, "-V", stdout = TRUE, stderr = TRUE, timeout = 20),
            error = function(e) character(0)))
        if (any(grepl("version", out, ignore.case = TRUE))) {
            ## A binary is there and runs. Is it the one this package expects?
            ## The build stamps DESCRIPTION's Version into the binary, and R
            ## normalises that same string for packageVersion(), so the two
            ## are directly comparable: "26.9.3" from `inla -V` must equal
            ## packageVersion("INLA"). They drift when a user upgrades the R
            ## package and keeps an older cached binary, which then fails in
            ## ways that look like modelling errors rather than a version
            ## mismatch. Report it; do not act on it.
            bv <- sub(".*version:[[:space:]]*", "",
                      grep("version", out, ignore.case = TRUE, value = TRUE)[1])
            bv <- trimws(bv)
            ## Compare against the binary version this package DECLARES it
            ## needs (Config/INLA/BinaryVersion), not against its own Version.
            ## The R code moves independently: most edits here need no new
            ## solver, so requiring equal versions would cry wolf on every R
            ## update. Only an OLDER binary than the declared minimum is a
            ## problem; a newer one is fine and stays quiet.
            need <- tryCatch(utils::packageDescription("INLA")[["Config/INLA/BinaryVersion"]],
                             error = function(e) NULL)
            if (!is.null(need) && nzchar(need) && nzchar(bv)) {
                older <- tryCatch(package_version(bv) < package_version(need),
                                  error = function(e) FALSE)
                if (isTRUE(older)) {
                    packageStartupMessage(
                        " - Binary is ", bv, " but this package needs ", need,
                        " or newer; run inla.stiles.install() to update it.")
                }
            }
            return(invisible(NULL))
        }
    }

    if (isTRUE(getOption("inla.stiles.autoinstall")) ||
        nzchar(Sys.getenv("INLA_STILES_AUTOINSTALL"))) {
        ## Opt-in only. Failure is reported and the session continues: an
        ## unreachable network must not stop library(INLA) from loading.
        res <- tryCatch(inla.stiles.install(), error = function(e) e)
        if (inherits(res, "error")) {
            packageStartupMessage("Could not install a binary: ", conditionMessage(res))
            packageStartupMessage("Run inla.stiles.install() to retry.")
        }
        return(invisible(NULL))
    }

    packageStartupMessage(
        " - No inla binary is installed yet; run inla.stiles.install() to fetch one.")
    invisible(NULL)
}

.onAttach <- function(...) {
    ## R CMD INSTALL starts helper R sessions (byte-compile, package indices,
    ## the load test) and each of them runs the user's ~/.Rprofile. A profile
    ## containing library(INLA) therefore attaches INLA three times while the
    ## package's own files are being replaced underneath it. Reading options
    ## out of a half-written lazy-load DB is not something this code can make
    ## safe, and it has segfaulted in exactly that window. So stay out of it
    ## entirely: there is nothing useful to say or do during an install anyway.
    ##
    ## R_CMD is the discriminator, verified by dumping the environment from a
    ## .Rprofile in both cases: it is set for any "R CMD ..." and absent in a
    ## normal session. R_INSTALL_PKG is NOT set at profile-read time, so it
    ## cannot be used here (nor in a .Rprofile, where it is often suggested).
    if (nzchar(Sys.getenv("R_CMD"))) {
        return(invisible(NULL))
    }
    if (interactive()) {
        inla.print.version()
    } else {
        packageStartupMessage(appendLF=FALSE)
    }
    try(inla.first.run.binary(), silent = TRUE)
}

.onUnload <- function(libpath) {
    ## library.dynam.unload("INLA", libpath)
}
