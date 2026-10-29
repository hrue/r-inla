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

.onAttach <- function(...) {
    if (nzchar(Sys.getenv("R_CMD"))) {
        return(invisible(NULL))
    }
    if (interactive()) {
        inla.print.version()
        if (is.null(inla.binary.path(check = FALSE))) {
            packageStartupMessage(" - No inla binary is installed yet; run inla.stiles.install() to fetch one.")
        }
    } else {
        packageStartupMessage(appendLF=FALSE)
    }
}

.onUnload <- function(libpath) {
    ## library.dynam.unload("INLA", libpath)
}
