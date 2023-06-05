#' View documentation
#' 
#' View documentation of latent, prior and likelihood models.
#' 
#' 
#' @aliases inla.doc inla.doc inla.inla.doc
#' @param what What to view documentation about; name of latent model, name of
#' prior, etc. (A regular expression.)
#' @param section An optional section, like \code{names(inla.models())}, to
#' look for the documentation. If missing, all sections are used.
#' @param verbose Logical If \code{TRUE} then run in verbose mode
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso \code{www.r-inla.org}
#' @examples
#' 
#' \dontrun{inla.doc("rw2")}
#' \dontrun{inla.doc("gaussian", section = "prior")}
#' 
#' @rdname doc
#' @export inla.doc
`inla.doc` <- function(what, section, verbose = FALSE) {
    what <- inla.trim.family(what)
    sections <- names(inla.models())
    m <- inla.models()

    for (sec in 1:length(m)) {
        if (missing(section) || sections[sec] == section) {
            if (verbose) {
                cat(paste("Search in section:", sections[sec], "\n"))
                cat(paste("\t", "Names in section:", names(m[[sec]]), "\n"))
            }
            x <- inla.trim.family(names(m[[sec]]))
            if (length(x) > 0L) {
                idx <- grep(what, x)
                if (length(idx) > 0L) {
                    for (i in idx) {
                        pdf <- m[[sec]][[i]]$pdf
                        if (!is.null(pdf) && !is.na(pdf)) {
                            filename <- paste("documentation/", sections[[sec]], "/", pdf, ".pdf", sep = "")
                            filename.full <- system.file(filename, package = "INLA")
                            if (verbose) {
                                cat("\n")
                                cat(paste("\t\t", "Found name:", names(m[[sec]])[i], "\n"))
                                cat(paste("\t\t", " ", "pdf           : ", pdf, ".pdf", sep = "", "\n"))
                                cat(paste("\t\t", " ", "filename      : ", filename, sep = "", "\n"))
                                cat(paste("\t\t", " ", "filename.full : ", filename.full, sep = "", "\n"))
                            }
                            if (file.exists(filename.full)) {
                                RShowDoc(
                                    paste("documentation/", sections[[sec]], "/", pdf, sep = ""),
                                    "pdf", "INLA"
                                )
                            }
                        }
                    }
                }
            }
        } else {
            if (verbose) {
                cat(paste("Skip section:", sections[sec], "\n"))
            }
        }
    }
}
