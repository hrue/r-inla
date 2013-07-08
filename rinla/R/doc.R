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
##!            name of latent model,  name of prior,  etc.}
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{\code{http://www.r-inla.org}}
##!\examples{
##!inla.doc("rw2")
##!inla.doc("gaussian")
##!}

`inla.doc` = function(what, verbose=FALSE)
{
    sections = names(inla.models())
    m = inla.models()
    for(section in 1:length(m)) {

        if (verbose) {
            cat(paste("Search in section:", names(m)[section], "\n"))
            cat(paste("\t",  "Names in section:",  names(m[[section]]), "\n"))
        }

        if (inla.is.element(what, m[[section]])) {
            idx = which(names(m[[section]]) == what)
            pdf = m[[section]][[idx]]$pdf

            if (verbose) {
                cat(paste("\t\t", "Found name:", names(m[[section]])[idx], "\n"))
                cat(paste("\t\t"," ", "Open pdf  : ", pdf, ".pdf", sep="", "\n"))
            }

            if (!is.null(pdf) && !is.na(pdf)) {
                ## just in case its more than one file
                for (pf in pdf) {
                    RShowDoc(paste(sections[[section]], "/", pf,  sep=""),
                             "pdf", "INLA")
                }
            }
        }
    }
}
    

    
