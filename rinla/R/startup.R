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

        hello <- paste0(
            "This is INLA_", version,
            " built ", built, ".", "\n",
            " - See www.r-inla.org/contact-us for how to get help.\n", 
            " - List available models/likelihoods/etc with inla.list.models()\n", 
            " - Use inla.doc(<NAME>) to access documentation"
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

## First run after installing: the package ships launcher scripts, not a
## solver, so `inla()` cannot run anything until a binary is present. Offer to
## fetch one instead of letting the first model fail with a message about a
## missing program.
##
## It ASKS rather than downloading on its own. A bundle is 26-130 MB, and a
## library() call silently pulling that over a metered or offline connection
## is not something a user can undo. Set
##   options(inla.stiles.autoinstall = TRUE)   (or INLA_STILES_AUTOINSTALL=1)
## to skip the question, which is what an unattended or container setup wants.
##
## Asked once: the answer, either way, is remembered in the cache directory,
## so declining does not turn into a prompt on every session.
`inla.first.run.binary` <- function() {
    cache <- tools::R_user_dir("INLA", "cache")
    stamp <- file.path(cache, "binary-offer-made")
    if (file.exists(stamp)) return(invisible(NULL))

    ## Already usable? Then say nothing: this is only for a fresh install.
    call <- tryCatch(inla.getOption("inla.call"), error = function(e) NULL)
    if (!is.null(call) && is.character(call) && nzchar(call) && file.exists(call)) {
        return(invisible(NULL))
    }

    auto <- isTRUE(getOption("inla.stiles.autoinstall")) ||
            nzchar(Sys.getenv("INLA_STILES_AUTOINSTALL"))

    if (!auto && !interactive()) {
        packageStartupMessage(
            "No inla binary is installed yet. Run inla.stiles.install() to fetch one.")
        return(invisible(NULL))
    }

    if (!auto) {
        packageStartupMessage("No inla binary is installed yet.")
        ans <- tryCatch(
            readline("Download and install one now? [y/N] "),
            error = function(e) "")
        ## Record the offer BEFORE acting: a failed or refused install must not
        ## leave the question to be asked again at every startup.
        dir.create(cache, recursive = TRUE, showWarnings = FALSE)
        try(writeLines(format(Sys.time()), stamp), silent = TRUE)
        if (!grepl("^[Yy]", ans)) {
            packageStartupMessage("Skipped. Run inla.stiles.install() when you want it.")
            return(invisible(NULL))
        }
    } else {
        dir.create(cache, recursive = TRUE, showWarnings = FALSE)
        try(writeLines(format(Sys.time()), stamp), silent = TRUE)
    }

    ## Never let this break library(INLA): a download failure is reported and
    ## the session continues, exactly as it would have without the offer.
    res <- tryCatch(inla.stiles.install(), error = function(e) e)
    if (inherits(res, "error")) {
        packageStartupMessage("Could not install a binary: ", conditionMessage(res))
        packageStartupMessage("Run inla.stiles.install() to retry.")
    }
    invisible(NULL)
}

.onAttach <- function(...) {
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
