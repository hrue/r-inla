## Export: inla.show.hyperspec

##!\name{inla.show.hyperspec}
##!\alias{inla.show.hyperspec}
##!\alias{show.hyperspec}
##!
##!\title{List all priors for the hyperparameters used}
##!
##!\description{List all priors for the hyperparameters used in a an \code{inla}-object}
##!
##!\usage{
##!inla.show.hyperspec(result)
##!}
##!\arguments{
## \item{result}{An \code {inla} -object} 
##!}
##!
##!\value{%%
##! Nothing. The used priors for the hyperparameters are displayed only.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\examples{
##!r = inla(y~1, data = data.frame(y=1:10))
##!inla.show.hyperspec(r)
##!}

`inla.show.hyperspec` = function(result)
{
    stopifnot(any(inherits(result, "inla")))
    tfile = tempfile()
    capture.output(str(result$all.hyper), file=tfile)
    all.hyper = readLines(tfile)
    unlink(tfile)
    
    all.hyper = gsub("\\.\\.", "  ", all.hyper)    
    for(r in c("inla\\.read\\.only",
               "attr\\(", "to\\.theta",
               "from\\.theta")) {
        idx = grep(r, all.hyper)
        if (length(idx) > 0) {
            all.hyper = all.hyper[-idx]
        }
    }

    cat(all.hyper, sep="\n")
    return (invisible())
}
