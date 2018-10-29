## Export: inla.debug.graph

##!\name{debug.graph}
##!\alias{debug.graph}
##!\alias{inla.debug.graph}
##!\title{Debug a graph-file}
##!\description{Debug a graph specification on file (ascii-mode only), by checking the specification along the way.}
##!\usage{inla.debug.graph(graph.file)}
##!\arguments{
##!    \item{graph.file}{The filename of the graph (ascii-mode)}
##!}
##!\value{
##! If an error is found, then an error message is shows, otherwise the graph-object returned by
##! \code{inla.read.graph()} is returned.
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\seealso{inla.read.graph}
##!\examples{
##!## Not run:
##!## cat("3\n 1 1 2n\ 2 1 1\n 3 4\n", file="g.dat")
##!## g = inla.debug.graph("g.dat")
##!## End(Not run)
##!}

`inla.debug.graph` = function(graph.file)
{
    ## read a graph with verbose output and try to detect any errors
    ## in the specification along the way. This is ment as a tool to
    ## detect errors in the graph specification only.

    stopifnot(file.exists(graph.file))

    xx = readLines(graph.file, encoding = "utf-8")
    
    cat("\n")
    cat("* File [", graph.file, "] consists of ", length(xx), " lines.\n", sep="")

    ## remove lines starting with '#'
    for(i in 1:length(xx)) {
        xx[i] = gsub("[ \t]+", " ", xx[i])
        xx[i] = gsub("#.*$", "", xx[i]) #
        if (length(grep("^[ \t]*$", xx[i])) > 0) {
            xx[i] = NA
        }
    }
    xx = xx[!is.na(xx)]
    cat("* Number of lines left after removing empty lines:", length(xx), "\n")

    to.ints = function(text) {
        xx = as.integer(unlist(sapply(text, function(x) strsplit(x, " "))))
        return (xx[!is.na(xx)])
    }
    
    N = to.ints(xx[1])
    stopifnot(N > 0)

    cat("* Size of the graph is N=", N, "\n")

    for(i in 2:length(xx)) {
        cat("* Read line", i, "...")

        x = to.ints(xx[i])
        cat("node =", x[1], " number.of.neigbours=", x[2])
        check = (x[2] + 2L == length(x))
        if (!check) {
            cat("\n*** ERROR IN THIS LINE: x = ", x, "\n")
            stop("Number of neigbours  (x[2]) does not correspond to the number of elements in that line")
        }
        
        check = all(sapply(x[-2], function(ii, N) return (ii >= 1L & ii <= N), N=N))
        if (!check) {
            cat("\n*** ERROR IN THIS LINE: x = ", x, "\n")
            if (min(x[-2]) == 0L) {
                stop("Minimum node is 0, which is not/no longer, allowed.")
            } else {
                stop("One of more of the nodes are (x[-2]) outside the legal range: [1, ...,", N, "]\n")
            }
        }
        if (x[2L] > 0L) {
            check = all(x[1L] != x[-c(1L, 2L)])
            if (!check) {
                cat("\n*** ERROR IN THIS LINE: x = ", x, "\n")
                stop("Node is defined as neighbour to itsef!")
            }
        }
        if (x[2L] > 0L) {
            check = (length(x) == length(unique(x[-c(1L, 2L)])) + 2L)
            if (!check) {
                cat("\n*** ERROR IN THIS LINE: x = ", x, "\n")
                stop("Some of the neighbours are duplicated!")
            }
        }
        cat("  ok\n")
    }

    cat("\n")
    cat("\n")
    cat("* I will now try to read the graph properly using inla.read.graph().\n")
    cat("* If there are any errors in the following, then recall\n")
    cat("* that the numbering for lines and nodes, are 0-based (and NOT 1-based)!\n")

    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.call.builtin()), "-s -m graph", shQuote(graph.file)), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.call.builtin()), "-s -m graph", shQuote(graph.file)), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    return  (inla.read.graph(graph.file))
}
