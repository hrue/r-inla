## Export: inla.list.models

##!\name{inla.list.models}
##!\alias{inla.list.models}
##!\alias{list.models}
##!\title{List available model components, likelihoods,  priors,  etc}
##!\description{
##!    List available model components, likelihoods,  priors,  etc. 
##!    To read specific documentation for the individual elements, use \code{\link{inla.doc}}.
##!}
##!\usage{
##!    inla.list.models(section = names(inla.models()), ...)
##!}
##!\arguments{
##!   \item{section}{The section(s) to list, missing \code{section} will list all sections.
##!                  \code{names(inla.models())} lists available sections.}
##!   \item{...}{Additional argument to \code{cat}}
##!}
##!\details{
##!   The list is \code{cat}'ed with \code{...} arguments.
##!
##!   This function is EXPERIMENTAL.
##!}
##!\value{
##!   Nothing is returned
##!}
##!\author{Havard Rue} 
##!\examples{
##!\dontrun{
##!inla.list.models("likelihood")
##!inla.list.models(c("prior", "group"))
##!inla.list.models(file=file("everything.txt"))
##!
##!#Show detailed doc for a specific prior/likelihood/latent model
##!inla.doc("binomial")
##!}
##!}


`inla.list.models` = function(section = names(inla.models()), ...)
{
    section = sort(match.arg(section, several.ok = TRUE))
    res.tot = NULL
    prefix = "    "
    for(sec in section) {
        res.tot = c(res.tot, paste0("Section [", sec, "]\n"))
        res = NULL
        a = inla.models()[[sec]]
        nm = names(a)
        for(i in seq_along(a)) {
            res = c(res, sprintf("%s%-30s%-40s\n", prefix, nm[i], a[[i]]$doc))
        }
        res.tot = c(res.tot, res[order(nm)])
    }
    cat(res.tot, ...)
    return(invisible())
}
