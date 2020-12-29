## Export: inla.version

## !\name{inla.version}
## !\alias{version}
## !\alias{inla.version}
## !
## !\title{Show the version of the INLA-package}
## !
## !\description{Show the version of the INLA-package}
## !
## !\usage{
## !inla.version(what = c("default", "version", "date"))
## !}
## !
## !\arguments{
## !  \item{what}{What to show version of}
## !}
## !
## !\value{%%
## !  \code{inla.version} display the current version information using \code{cat}
## !  with
## !  \code{default} or \code{info},  or return
## !  other spesific requests through the call.
## !}
## !%%
## !
## !\author{Havard Rue \email{hrue@r-inla.org}}
## !
## !\examples{
## !## Summary of all
## !inla.version()
## !## The building date
## !inla.version("bdate")
## !}

`inla.version` <- function(what = c("default", "version", "date")) {
    `trim` <- function(string) {
        string <- gsub("^[ \t]+", "", string)
        string <- gsub("[ \t]+$", "", string)
        return(string)
    }

    date <- trim("Thu 30 Jun 08:30:28 CEST 2016")
    version <- trim("Version_12.34.56")
    what <- match.arg(what)

    if (what %in% "default") {
        cat("\n")
        cat(paste("\n\tR-INLA version ..........: ", version, "\n", sep = ""))
        cat(paste("\tDate ....................: ", date, "\n", sep = ""))
        cat("\tMaintainers .............: Havard Rue <hrue@r-inla.org>\n")
        cat("\t                         : Finn Lindgren <finn.lindgren@gmail.com>\n")
        cat("\t                         : Elias Teixeira Krainski <elias.krainski@gmail.com>\n")
        cat("\tMain web-page ...........: www.r-inla.org\n")
        cat("\tDownload-page ...........: inla.r-inla-download.org\n")
        cat("\tEmail support ...........: help@r-inla.org\n")
        cat("\t                         : r-inla-discussion-group@googlegroups.com\n")
        cat("\tSource-code .............: github.com/hrue/r-inla\n")
        cat("\n")

        return(invisible())
    } else if (what %in% "date") {
        return(date)
    } else if (what %in% "version") {
        return(version)
    }

    stop("This should not happen.")
}
