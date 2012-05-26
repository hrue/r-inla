##! \name{qreordering}
##! \alias{inla.qreordering}
##! \alias{qreordering}
##! 
##! \title{Compute the reordering using the GMRFLib implementation}
##! 
##! \description{This function compute the reordering (or find the best reordering)
##!              using the GMRFLib implementation}
##! \usage{
##!     inla.qreordering(graph, reordering)
##! }
##! 
##! \arguments{
##!   \item{graph}{A (inla-)graph object, a filename containing the graph or a matrix defining it.}
##!   \item{reordering}{The type of reordering algorithm,  one of
##!        "auto", "default", "identity", "band", "metis", "genmmd", "amd", "amdbar", "md", "mmd", "amdc" and "amdbarc".
##!        The default is "auto" which try several and find the best one.}
##!}
##!\value{
##!  \code{inla.qreordering} returns a list with the name of the reordering algorithm used or found,  the reordering code for the reordering algorithm,
##!                          the actual reordering and its inverse.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##! g = system.file("demodata/germany.graph", package="INLA")
##! r = inla.qreordering(g)
##!}

`inla.qreordering` = function(graph,
        reordering = c("auto", "default", "identity", "band", "metis", "genmmd", "amd", "amdbar", "md", "mmd", "amdc", "amdbarc"))
{
    reordering = match.arg(reordering)

    g = inla.read.graph(graph)
    g.file = inla.write.graph(g)

    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qreordering", 
                "-r", reordering, g.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qreordering",
                "-r", reordering, g.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    unlink(g.file)

    nm = as.character(s[1L])
    code = as.numeric(s[2L])
    s = as.numeric(s[-c(1L, 2L)])
    r = s + 1L
    ir = r
    ir[r] = 1L:length(r)

    return ( list(name = nm, code = code, reordering = r, ireordering = ir) )
}
