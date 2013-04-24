`inla.version` = function (what = c("default",
                                   "info", 
                                   "hgid",
                                   "rinla",
                                   "inla",
                                   "date",
                                   "bdate"), 
        ## for backward compatibility
        details = FALSE, quiet = FALSE, hgid = FALSE) 
{
    rinla.hgid = "EXAMPLE hgid: bc0bd174e33d  date: Mon Mar 05 10:05:06 2012 +0100"
    inla.hgid = "EXAMPLE hgid: bc0bd174e33d  date: Mon Mar 05 10:05:06 2012 +0100"
    date = "EXAMPLE Mon Mar 5 10:15:43 CET 2012"
    bdate = "EXAMPLE201203071055"
    
    if (missing(what) && (!missing(details) || !missing(quiet) || !missing(hgid))) {
        ## these options were added for backward compatibility: will
        ## not stay forever.
        if (hgid) {
            return (inla.version("hgid"))
        }
        if (!quiet) {
            return (inla.version("default"))
        } else {
            return (inla.version("bdate"))
        }
        stop("This should not happen.")
    }

    what = match.arg(what)

    if (what %in% c("default", "info")) {

        cat("\n")
        cat(paste("\n\tINLA build date .........: ",  date, "\n",  sep=""))
        cat(paste(  "\tINLA hgid ...............: ", rinla.hgid, "\n", sep=""))
        cat(paste(  "\tINLA-program hgid .......: ", inla.hgid, "\n", sep=""))
        cat(        "\tMaintainers .............: Havard Rue <hrue@math.ntnu.no>\n")
        cat(        "\t                         : Finn Lindgren <finn.lindgren@gmail.com>\n")
        cat(        "\t                         : Daniel Simpson <dp.simpson@gmail.com>\n")
        cat(        "\t                         : Andrea Riebler <andrea.riebler@math.ntnu.no>\n")
        cat(        "\tWeb-page ................: http://www.r-inla.org\n")
        cat(        "\tEmail support ...........: help@r-inla.org\n")
        cat(        "\t                         : r-inla-discussion-group@googlegroups.com\n")
        cat(        "\tSource-code .............: http://inla.googlecode.com\n")
        cat("\n")

        return(invisible())

    } else if (what %in% c("hgid", "rinla")) {
        return (rinla.hgid)
    } else if (what %in% "inla") {
        return (inla.hgid)
    } else if (what %in% "date") {
        return (date) 
    } else if (what %in% "bdate") {
        return (bdate) 
    }
 
    stop("This should not happen.")
}
