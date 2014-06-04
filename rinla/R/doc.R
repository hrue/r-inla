## Export: inla.doc

##!\name{inla.doc}
##!\alias{inla.doc}
##!\alias{inla.doc}
##!\alias{inla.inla.doc}
##!\title{View documentation}
##!\description{View documentation of latent, prior and likelihood models.}
##!\usage{
##!inla.doc(what, verbose=FALSE)
##!}
##!\arguments{
##!\item{what}{What to view documentation about;
##!            name of latent model,  name of prior,  etc. (A regular expression.)}
##!\item{verbose}{Logical if \code{TRUE} then run in verbose mode}
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{\code{http://www.r-inla.org}}
##!\examples{
##!\dontrun{inla.doc("rw2")}
##!\dontrun{inla.doc("gaussian")}
##!}

`inla.doc` = function(what, verbose=FALSE)
{
    what = inla.trim.family(what)
    sections = names(inla.models())
    m = inla.models()

    for(section in 1:length(m)) {
        if (verbose) {
            cat(paste("Search in section:", names(m)[section], "\n"))
            cat(paste("\t",  "Names in section:",  names(m[[section]]), "\n"))
        }

        x = inla.trim.family(names(m[[section]]))
        if (length(x) > 0L) {
            idx = grep(what, inla.trim.family(names(m[[section]])))
        } else {
            idx = numeric(0)
        }
        if (length(idx) > 0L) {
            for(i in idx) {
                pdf = m[[section]][[i]]$pdf
                if (!is.null(pdf) && !is.na(pdf)) {
                    if (verbose) {
                        cat(paste("\t\t", "Found name:", names(m[[section]])[i], "\n"))
                        cat(paste("\t\t"," ", "Open pdf  : ", pdf, ".pdf", sep="", "\n"))
                    }
                    RShowDoc(paste(sections[[section]], "/", pdf,  sep=""),
                             "pdf", "INLA")
                }
            }
        }
    }
}



