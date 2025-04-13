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

        if (FALSE) {
            if (!inla.os("windows")) {
                hello <- paste0(hello,
                                "\n",
                                " - To enable PARDISO sparse library; see inla.pardiso()"
                                )
            }
        }
        
        s <- inla.prune.check()
        if (s > 0) {
            hello <- paste0(
                hello,
                "\n",
                paste0(
                    " - Save ", s,
                    "Mb of storage running 'inla.prune()'"
                )
            )
        }

        opts <- options()
        options(timeout = 2)
        vers <- try(readLines("https://inla.r-inla-download.org/VERSIONS",
                              n = 4, encoding = "UTF-8"), silent = TRUE)
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

        packageStartupMessage(hello)
    }
}

.onLoad <- function(...) {
    ## nothing for the moment
}

.onAttach <- function(...) {
    inla.print.version()
}

.onUnload <- function(libpath) {
    ## library.dynam.unload("INLA", libpath)
}
