## Export: inla.doc

##!\name{inla.doc}
##!\alias{inla.doc}
##!\alias{inla.doc}
##!\alias{inla.inla.doc}
##!\title{View documentation}
##!\description{View documentation of latent, prior and likelihood models.}
##!\usage{
##!inla.doc(what, sec, verbose=FALSE)
##!}
##!\arguments{
##!\item{what}{What to view documentation about;
##!            name of latent model,  name of prior,  etc. (A regular expression.)}
##!\item{sec}{An optional section to look for the documentation. If missing, all sections are used.}
##!\item{verbose}{Logical if \code{TRUE} then run in verbose mode}
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\seealso{\code{www.r-inla.org}}
##!\examples{
##!\dontrun{inla.doc("rw2")}
##!\dontrun{inla.doc("gaussian")}
##!}

`inla.doc` = function(what, sec, verbose=FALSE)
{
    what = inla.trim.family(what)
    sections = names(inla.models())
    m = inla.models()

    for(section in 1:length(m)) {
        if (missing(sec) || sections[section] == sec) {
            if (verbose) {
                cat(paste("Search in section:", sections[section], "\n"))
                cat(paste("\t",  "Names in section:",  names(m[[section]]), "\n"))
            }
            x = inla.trim.family(names(m[[section]]))
            if (length(x) > 0L) {
                idx = grep(what, x)
                if (length(idx) > 0L) {
                    for(i in idx) {
                        pdf = m[[section]][[i]]$pdf
                        if (!is.null(pdf) && !is.na(pdf)) {
                            filename = paste("doc/", sections[[section]], "/", pdf, ".pdf", sep="")
                            filename.full = system.file(filename, package="INLA")
                            if (verbose) {
                                cat("\n")
                                cat(paste("\t\t", "Found name:", names(m[[section]])[i], "\n"))
                                cat(paste("\t\t"," ", "pdf           : ", pdf, ".pdf", sep="", "\n"))
                                cat(paste("\t\t"," ", "filename      : ", filename, sep="", "\n"))
                                cat(paste("\t\t"," ", "filename.full : ", filename.full, sep="", "\n"))
                            }
                            if (file.exists(filename.full)) {
                                RShowDoc(paste(sections[[section]], "/", pdf,  sep=""),
                                         "pdf", "INLA")
                            }
                        }
                    }
                }
            }
        } else {
            if (verbose) {
                cat(paste("Skip section:", sections[section], "\n"))
            }
        }
    }
}





