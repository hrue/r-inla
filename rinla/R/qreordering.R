#' Compute the reordering using the GMRFLib implementation
#' 
#' This function compute the reordering (or find the best reordering) using the
#' GMRFLib implementation
#' 
#' 
#' @aliases inla.qreordering qreordering
#' @param graph A \code{(inla-)graph} object
#' @param reordering The name of the reordering algorithm to be used; either
#' one of the names listed in \code{inla.reorderings()}.  The default is "auto"
#' which try several reordering algorithm and use the best one for this
#' particular matrix.
#' @return \code{inla.qreordering} returns a list with the name of the
#' reordering algorithm used or found, the reordering code for the reordering
#' algorithm, the actual reordering and its inverse.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'  g = system.file("demodata/germany.graph", package="INLA")
#'  r = inla.qreordering(g)
#'  m = inla.graph2matrix(g)
#'  r = inla.qreordering(m)
#' 
#' @name qreordering
#' @rdname qreordering
#' @export

`inla.qreordering` <- function(graph, reordering = inla.reorderings()) {
    reordering <- match.arg(reordering)

    if (inla.is.matrix(graph)) {
        graph <- inla.sparse.check(graph)
        g.file <- inla.write.fmesher.file(graph)
        g.remove <- TRUE
    } else if (inla.is.fmesher.file(graph)) {
        g.file <- graph
        g.remove <- FALSE
    } else {
        g <- inla.read.graph(graph)
        g.file <- inla.write.graph(g)
        g.remove <- FALSE
    }

    ## smtp must be taucs
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qreordering",
            "-r", reordering, "-S", "taucs", g.file
        ), intern = TRUE)
    } else if (inla.os("windows")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qreordering",
            "-r", reordering, "-S", "taucs", g.file
        ), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    if (g.remove) {
        unlink(g.file)
    }

    for (i in 1:length(s)) {
        ## the output containts this code, as there may be some other output before this
        if (as.character(s[i]) == "QREORDERING") {
              break
          }
    }
    s <- s[-(1:i)]

    nm <- as.character(s[1L])
    code <- as.numeric(s[2L])
    s <- as.numeric(s[-c(1L, 2L)])
    r <- s + 1L
    ir <- r
    ir[r] <- 1L:length(r)

    return(list(name = nm, code = code, reordering = r, ireordering = ir))
}
